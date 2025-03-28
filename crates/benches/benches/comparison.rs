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
    generate_comparison!(11, c);
    generate_comparison!(12, c);
    generate_comparison!(13, c);
    generate_comparison!(14, c);
    generate_comparison!(15, c);
    generate_comparison!(16, c);
    generate_comparison!(17, c);
    generate_comparison!(18, c);
    generate_comparison!(19, c);
    generate_comparison!(20, c);
    generate_comparison!(21, c);
    generate_comparison!(22, c);
    generate_comparison!(23, c);
    generate_comparison!(24, c);
    generate_comparison!(25, c);
    generate_comparison!(26, c);
    generate_comparison!(27, c);
    generate_comparison!(28, c);
    generate_comparison!(29, c);
    generate_comparison!(30, c);
    generate_comparison!(31, c);
    generate_comparison!(32, c);
    generate_comparison!(33, c);
    generate_comparison!(34, c);
    generate_comparison!(35, c);
    generate_comparison!(36, c);
    generate_comparison!(37, c);
    generate_comparison!(38, c);
    generate_comparison!(39, c);
    generate_comparison!(40, c);
    generate_comparison!(41, c);
    generate_comparison!(42, c);
    generate_comparison!(43, c);
    generate_comparison!(44, c);
    generate_comparison!(45, c);
    generate_comparison!(46, c);
    generate_comparison!(47, c);
    generate_comparison!(48, c);
    generate_comparison!(49, c);
    generate_comparison!(50, c);
    generate_comparison!(51, c);
    generate_comparison!(52, c);
    generate_comparison!(53, c);
    generate_comparison!(54, c);
    generate_comparison!(55, c);
    generate_comparison!(56, c);
    generate_comparison!(57, c);
    generate_comparison!(58, c);
    generate_comparison!(59, c);
    generate_comparison!(60, c);
    generate_comparison!(61, c);
    generate_comparison!(62, c);
    generate_comparison!(63, c);
    generate_comparison!(64, c);
    generate_comparison!(65, c);
    generate_comparison!(66, c);
    generate_comparison!(67, c);
    generate_comparison!(68, c);
    generate_comparison!(69, c);
    generate_comparison!(70, c);
    generate_comparison!(71, c);
    generate_comparison!(72, c);
    generate_comparison!(73, c);
    generate_comparison!(74, c);
    generate_comparison!(75, c);
    generate_comparison!(76, c);
    generate_comparison!(77, c);
    generate_comparison!(78, c);
    generate_comparison!(79, c);
    generate_comparison!(80, c);
    generate_comparison!(81, c);
    generate_comparison!(82, c);
    generate_comparison!(83, c);
    generate_comparison!(84, c);
    generate_comparison!(85, c);
    generate_comparison!(86, c);
    generate_comparison!(87, c);
    generate_comparison!(88, c);
    generate_comparison!(89, c);
    generate_comparison!(90, c);
    generate_comparison!(91, c);
    generate_comparison!(92, c);
    generate_comparison!(93, c);
    generate_comparison!(94, c);
    generate_comparison!(95, c);
    generate_comparison!(96, c);
    generate_comparison!(97, c);
    generate_comparison!(98, c);
    generate_comparison!(99, c);
    generate_comparison!(100, c);
    generate_comparison!(101, c);
    generate_comparison!(102, c);
    generate_comparison!(103, c);
    generate_comparison!(104, c);
    generate_comparison!(105, c);
    generate_comparison!(106, c);
    generate_comparison!(107, c);
    generate_comparison!(108, c);
    generate_comparison!(109, c);
    generate_comparison!(110, c);
    generate_comparison!(111, c);
    generate_comparison!(112, c);
    generate_comparison!(113, c);
    generate_comparison!(114, c);
    generate_comparison!(115, c);
    generate_comparison!(116, c);
    generate_comparison!(117, c);
    generate_comparison!(118, c);
    generate_comparison!(119, c);
    generate_comparison!(120, c);
    generate_comparison!(121, c);
    generate_comparison!(122, c);
    generate_comparison!(123, c);
    generate_comparison!(124, c);
    generate_comparison!(125, c);
    generate_comparison!(126, c);
    generate_comparison!(127, c);
    generate_comparison!(128, c);
    generate_comparison!(129, c);
    generate_comparison!(130, c);
    generate_comparison!(131, c);
    generate_comparison!(132, c);
    generate_comparison!(133, c);
    generate_comparison!(134, c);
    generate_comparison!(135, c);
    generate_comparison!(136, c);
    generate_comparison!(137, c);
    generate_comparison!(138, c);
    generate_comparison!(139, c);
    generate_comparison!(140, c);
    generate_comparison!(141, c);
    generate_comparison!(142, c);
    generate_comparison!(143, c);
    generate_comparison!(144, c);
    generate_comparison!(145, c);
    generate_comparison!(146, c);
    generate_comparison!(147, c);
    generate_comparison!(148, c);
    generate_comparison!(149, c);
    generate_comparison!(150, c);
    generate_comparison!(151, c);
    generate_comparison!(152, c);
    generate_comparison!(153, c);
    generate_comparison!(154, c);
    generate_comparison!(155, c);
    generate_comparison!(156, c);
    generate_comparison!(157, c);
    generate_comparison!(158, c);
    generate_comparison!(159, c);
    generate_comparison!(160, c);
    generate_comparison!(161, c);
    generate_comparison!(162, c);
    generate_comparison!(163, c);
    generate_comparison!(164, c);
    generate_comparison!(165, c);
    generate_comparison!(166, c);
    generate_comparison!(167, c);
    generate_comparison!(168, c);
    generate_comparison!(169, c);
    generate_comparison!(170, c);
    generate_comparison!(171, c);
    generate_comparison!(172, c);
    generate_comparison!(173, c);
    generate_comparison!(174, c);
    generate_comparison!(175, c);
    generate_comparison!(176, c);
    generate_comparison!(177, c);
    generate_comparison!(178, c);
    generate_comparison!(179, c);
    generate_comparison!(180, c);
    generate_comparison!(181, c);
    generate_comparison!(182, c);
    generate_comparison!(183, c);
    generate_comparison!(184, c);
    generate_comparison!(185, c);
    generate_comparison!(186, c);
    generate_comparison!(187, c);
    generate_comparison!(188, c);
    generate_comparison!(189, c);
    generate_comparison!(190, c);
    generate_comparison!(191, c);
    generate_comparison!(192, c);
    generate_comparison!(193, c);
    generate_comparison!(194, c);
    generate_comparison!(195, c);
    generate_comparison!(196, c);
    generate_comparison!(197, c);
    generate_comparison!(198, c);
    generate_comparison!(199, c);
    generate_comparison!(200, c);
}

criterion_group!(benches, compare);
criterion_main!(benches);
