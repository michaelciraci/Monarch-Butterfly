use monarch_butterfly::{fft, ifft};
use num_complex::Complex;
use num_traits::{Float, FloatConst};
use rustfft::{FftNum, FftPlanner};

fn check<const N: usize, T: Float + FloatConst + FftNum>(tolerance: T) {
    let dense: [Complex<T>; N] = core::array::from_fn(|i| {
        Complex::new(
            T::from_f64((i * 13 % 17) as f64 / 8.0 - 1.0).unwrap(),
            T::from_f64((i * i % 11) as f64 / 7.0 - 0.5).unwrap(),
        )
    });
    let zero = Complex::new(T::zero(), T::zero());
    let mut planner = FftPlanner::<T>::new();
    let forward = planner.plan_fft_forward(N);
    let inverse = planner.plan_fft_inverse(N);

    // An impulse at every position checks indexing and natural bin order.
    for impulse_index in 0..=N {
        let input = if impulse_index == N {
            dense
        } else {
            let mut impulse = [zero; N];
            impulse[impulse_index] = Complex::new(T::one(), -T::one());
            impulse
        };
        for (actual, plan) in [
            (fft::<N, _, _>(input), &forward),
            (ifft::<N, _, _>(&input[..]), &inverse),
        ] {
            let mut expected = input;
            plan.process(&mut expected);
            for (k, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
                assert!(
                    (actual - expected).norm() < tolerance,
                    "size {N}, bin {k}, direction {:?}: {actual:?} != {expected:?}",
                    plan.fft_direction()
                );
            }
        }
    }
}

#[test]
fn prime_kernels_and_composition_match_rustfft_for_varied_complex_inputs() {
    fn check_many<const N: usize, T: Float + FloatConst + FftNum>(tolerance: T) {
        let mut planner = FftPlanner::<T>::new();
        let forward = planner.plan_fft_forward(N);
        let inverse = planner.plan_fft_inverse(N);
        let mut state = 0x9e37_79b9_u32;
        for _ in 0..128 {
            let mut sample = || {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                T::from_f64(f64::from(state) / f64::from(u32::MAX) * 2.0 - 1.0).unwrap()
            };
            let input: [Complex<T>; N] = core::array::from_fn(|_| Complex::new(sample(), sample()));
            for (actual, plan) in [
                (fft::<N, _, _>(input), &forward),
                (ifft::<N, _, _>(input), &inverse),
            ] {
                let mut expected = input;
                plan.process(&mut expected);
                for (bin, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
                    assert!(
                        (actual - expected).norm() < tolerance,
                        "size {N}, bin {bin}, direction {:?}: {actual:?} != {expected:?}",
                        plan.fft_direction()
                    );
                }
            }
        }
    }

    check_many::<5, f32>(1e-5);
    check_many::<5, f64>(1e-12);
    check_many::<7, f32>(1e-5);
    check_many::<7, f64>(1e-12);
    check_many::<10, f32>(1e-5);
    check_many::<10, f64>(1e-12);
}

macro_rules! check_size {
    ($name:ident, $n:literal) => {
        #[test]
        fn $name() {
            check::<$n, f32>(1e-5);
            check::<$n, f64>(1e-12);
        }
    };
}

check_size!(size_1, 1);
check_size!(size_2, 2);
check_size!(size_3, 3);
check_size!(size_4, 4);
check_size!(size_5, 5);
check_size!(size_6, 6);
check_size!(size_7, 7);
check_size!(size_8, 8);
check_size!(size_9, 9);
check_size!(size_10, 10);
