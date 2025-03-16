use butterfly::*;
use criterion::{Criterion, black_box, criterion_group, criterion_main};
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
            let _ = butterfly16(black_box(input));
        })
    });
}

fn compare32(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(32);

    let mut input = vec![Complex::ZERO; 32];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-32", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::ZERO; 32];

    c.bench_function("michael-32", |b| {
        b.iter(|| {
            let _ = butterfly32(black_box(input));
        })
    });
}

fn compare64(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(64);

    let mut input = vec![Complex::ZERO; 64];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-64", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::ZERO; 64];

    c.bench_function("michael-64", |b| {
        b.iter(|| {
            let _ = butterfly64(black_box(input));
        })
    });
}

fn compare128(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(128);

    let mut input = vec![Complex::ZERO; 128];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-128", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::ZERO; 64];

    c.bench_function("michael-128", |b| {
        b.iter(|| {
            let _ = butterfly64(black_box(input));
        })
    });
}

criterion_group!(benches, compare128);
criterion_main!(benches);
