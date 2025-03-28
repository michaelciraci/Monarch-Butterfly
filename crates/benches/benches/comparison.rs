use criterion::{Criterion, black_box, criterion_group, criterion_main};
use fftw::{
    array::AlignedVec,
    plan::{C2CPlan, C2CPlan32},
    types::{Flag, Sign},
};
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
                let _ = fft::<$idx, _, _>(black_box(input));
            })
        });

        let n = $idx;
        let mut plan: C2CPlan32 = C2CPlan::aligned(&[n], Sign::Forward, Flag::MEASURE).unwrap();
        let mut a = AlignedVec::new(n);
        let mut o = AlignedVec::new(n);

        $c.bench_function(&format!("fftw-{}", $idx), |b| {
            b.iter(|| {
                plan.c2c(&mut a, &mut o).unwrap();
            })
        });
    }};
}

fn compare(c: &mut Criterion) {
    generate_comparison!(9, c);
    generate_comparison!(10, c);
    generate_comparison!(11, c);
    generate_comparison!(12, c);
    generate_comparison!(13, c);
    generate_comparison!(14, c);
    generate_comparison!(15, c);
    generate_comparison!(16, c);
}

criterion_group!(benches, compare);
criterion_main!(benches);
