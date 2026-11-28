use criterion::{Criterion, black_box, criterion_group, criterion_main};
// use fftw::{
//     array::AlignedVec,
//     plan::{C2CPlan, C2CPlan32},
//     types::{Flag, Sign},
// };
use monarch_butterfly::*;
use num_complex::Complex;
use rustfft::FftPlanner;

macro_rules! generate_comparison {
    ($idx:expr, $c:ident) => {{
        let mut planner = FftPlanner::<f32>::new();
        let p = planner.plan_fft_forward($idx);

        let mut input = vec![Complex::<f32>::ZERO; $idx];
        let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

        $c.bench_function(&format!("rustfft-{}", $idx), |b| {
            b.iter(|| {
                p.process_with_scratch(&mut input, &mut scratch);
            })
        });

        let input = [Complex::<f32>::ZERO; $idx];

        $c.bench_function(&format!("monarch-{}", $idx), |b| {
            b.iter(|| {
                black_box(fft::<$idx, _, _>(black_box(input)));
            })
        });

        // let n = $idx;
        // let mut plan: C2CPlan32 = C2CPlan::aligned(&[n], Sign::Forward, Flag::MEASURE).unwrap();
        // let mut a = AlignedVec::new(n);
        // let mut o = AlignedVec::new(n);

        // $c.bench_function(&format!("fftw-{}", $idx), |b| {
        //     b.iter(|| {
        //         plan.c2c(&mut a, &mut o).unwrap();
        //     })
        // });
    }};
}

// Compare the existing by-value harness with borrowed inputs. Include dense
// nonzero samples because the compiler may lay out constant inputs differently.
fn input_diagnostics<const N: usize>(c: &mut Criterion) {
    let dense: [Complex<f32>; N] =
        core::array::from_fn(|i| Complex::new(i as f32 * 0.25 - 0.5, i as f32 * 0.1 - 0.3));
    for (label, input) in [("zero", [Complex::ZERO; N]), ("dense", dense)] {
        c.bench_function(&format!("monarch-diagnostic-{N}-{label}-value"), |b| {
            b.iter(|| black_box(fft::<N, _, _>(black_box(input))))
        });
        c.bench_function(&format!("monarch-diagnostic-{N}-{label}-borrowed"), |b| {
            b.iter(|| black_box(fft::<N, _, _>(black_box(&input))))
        });
    }
}

fn compare(c: &mut Criterion) {
    generate_comparison!(1, c);
    generate_comparison!(2, c);
    generate_comparison!(3, c);
    generate_comparison!(4, c);
    generate_comparison!(5, c);
    generate_comparison!(6, c);
    generate_comparison!(7, c);
    generate_comparison!(8, c);
    generate_comparison!(9, c);
    generate_comparison!(10, c);

    input_diagnostics::<5>(c);
    input_diagnostics::<10>(c);
}

criterion_group!(benches, compare);
criterion_main!(benches);
