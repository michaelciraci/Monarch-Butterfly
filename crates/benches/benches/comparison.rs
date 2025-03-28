use criterion::{Criterion, black_box, criterion_group, criterion_main};
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
    }};
}

fn compare(c: &mut Criterion) {
    generate_comparison!(79, c);
    generate_comparison!(80, c);
    generate_comparison!(81, c);
    generate_comparison!(82, c);
    generate_comparison!(83, c);
    generate_comparison!(84, c);
    generate_comparison!(85, c);
    generate_comparison!(86, c);
}

criterion_group!(benches, compare);
criterion_main!(benches);
