use criterion::{black_box, criterion_group, criterion_main, Criterion};
use monarch_butterfly::*;
use num_complex::Complex;
use rustfft::FftPlanner;

fn compare16(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(16);

    let mut input = vec![Complex::<f32>::ZERO; 16];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-16", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::<f32>::ZERO; 16];

    c.bench_function("monarch-16", |b| {
        b.iter(|| {
            let _ = fft16(black_box(input));
        })
    });
}

fn compare32(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(32);

    let mut input = vec![Complex::<f32>::ZERO; 32];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-32", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::<f32>::ZERO; 32];

    c.bench_function("monarch-32", |b| {
        b.iter(|| {
            let _ = fft32(black_box(input));
        })
    });
}

fn compare64(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(64);

    let mut input = vec![Complex::<f32>::ZERO; 64];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-64", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::<f32>::ZERO; 64];

    c.bench_function("monarch-64", |b| {
        b.iter(|| {
            let _ = fft64(black_box(input));
        })
    });
}

fn compare128(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(128);

    let mut input = vec![Complex::<f32>::ZERO; 128];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-128", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::<f32>::ZERO; 128];

    c.bench_function("monarch-128", |b| {
        b.iter(|| {
            let _ = fft128(black_box(input));
        })
    });
}

fn compare256(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(256);

    let mut input = vec![Complex::<f32>::ZERO; 256];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-256", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::<f32>::ZERO; 256];

    c.bench_function("monarch-256", |b| {
        b.iter(|| {
            let _ = fft256(black_box(input));
        })
    });
}

fn compare512(c: &mut Criterion) {
    let mut planner = FftPlanner::<f32>::new();
    let p = planner.plan_fft_forward(512);

    let mut input = vec![Complex::<f32>::ZERO; 512];
    let mut scratch = vec![Complex::ZERO; p.get_inplace_scratch_len()];

    c.bench_function("rustfft-512", |b| {
        b.iter(|| {
            p.process_with_scratch(&mut input, &mut scratch);
        })
    });

    let input = [Complex::<f32>::ZERO; 512];

    c.bench_function("monarch-512", |b| {
        b.iter(|| {
            let _ = fft512(black_box(input));
        })
    });
}

criterion_group!(benches, compare16, compare32, compare64, compare128, compare256, compare512);
criterion_main!(benches);
