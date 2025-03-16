use butterfly::butterfly16;
use criterion::{Criterion, criterion_group, criterion_main};
use num_complex::Complex;
use rustfft::FftPlanner;

fn compare16(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(16);

    let mut input = vec![Complex::ZERO; 16];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-16", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::ZERO; 16];

    c.bench_function("michael-16", |b| {
        b.iter(|| {
            let _ = butterfly16(input);
        })
    });
}

criterion_group!(benches, compare16);
criterion_main!(benches);
