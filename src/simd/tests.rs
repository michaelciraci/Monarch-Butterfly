extern crate std;

use core::{array::from_fn, fmt::Debug};
use num_traits::FloatConst;

use super::*;
use crate::{fft, fixed, ifft};

fn check<T: Float + FloatConst + Debug>(input: [Complex<T>; 8]) {
    let simd = [
        try_fft8::<8, false, _>(&input).unwrap(),
        try_fft8::<8, true, _>(&input).unwrap(),
    ];
    let scalar = [
        fixed::fft8::<false, _>(input),
        fixed::fft8::<true, _>(input),
    ];
    let public = [fft::<8, _, _>(&input), ifft::<8, _, _>(&input)];
    for ((simd, scalar), public) in simd.into_iter().zip(scalar).zip(public) {
        for ((a, b), c) in simd.into_iter().zip(scalar).zip(public) {
            for (a, b, c) in [(a.re, b.re, c.re), (a.im, b.im, c.im)] {
                assert!(
                    a == b || (a.is_nan() && b.is_nan()),
                    "SIMD {a:?}, scalar {b:?}"
                );
                assert!(
                    a == c || (a.is_nan() && c.is_nan()),
                    "SIMD {a:?}, public {c:?}"
                );
            }
        }
    }
}

fn check_inputs<T: Float + FloatConst + Debug>() {
    let zero = Complex::new(T::zero(), T::zero());
    for index in 0..8 {
        for value in [
            T::one(),
            -T::zero(),
            T::infinity(),
            T::neg_infinity(),
            T::nan(),
        ] {
            for impulse in [
                Complex::new(value, T::zero()),
                Complex::new(T::zero(), value),
            ] {
                let mut input = [zero; 8];
                input[index] = impulse;
                check(input);
            }
        }
    }

    let mut state = 0x619B_34F7_u32;
    let mut sample = || {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        T::from(f64::from(state) / 2_147_483_648.0 - 1.0).unwrap()
    };
    for scale in [
        T::min_positive_value(),
        T::one(),
        T::max_value() / T::from(64).unwrap(),
    ] {
        for _ in 0..128 {
            check(from_fn(|_| {
                Complex::new(sample() * scale, sample() * scale)
            }));
        }
    }
}

#[test]
fn both_precisions_match_scalar_for_impulses_special_values_and_random_inputs() {
    check_inputs::<f32>();
    check_inputs::<f64>();
}

fn unaligned<P: Pair<Scalar = T>, T: Float + Debug>() {
    #[repr(C, align(16))]
    struct Guarded<T>(T, [Complex<T>; 8], T);

    let guard = T::one();
    let input = Guarded(
        guard,
        from_fn(|i| Complex::new(T::from(i).unwrap(), -guard)),
        guard,
    );
    let mut output = Guarded(guard, [Complex::new(T::zero(), T::zero()); 8], guard);
    let src = ptr::from_ref(&input.1);
    let dst = ptr::from_mut(&mut output.1);
    assert_ne!(src.addr() % 16, 0);
    assert_ne!(dst.addr() % 16, 0);
    // SAFETY: The disjoint f32/f64 arrays each have sixteen scalar lanes.
    unsafe { write_fft8::<false, P>(src.cast(), dst.cast()) };
    assert_eq!(output.1, fixed::fft8::<false, _>(input.1));
    assert_eq!((output.0, output.2), (guard, guard));
}

#[test]
fn unaligned_loads_and_stores_preserve_guards() {
    unaligned::<arch::F32Pair, f32>();
    unaligned::<F64Pair, f64>();
}

#[test]
fn dispatch_validates_length() {
    let input = [Complex::new(0.0_f32, 0.0); 9];
    assert!(try_fft8::<9, false, _>(&input).is_none());
    assert!(std::panic::catch_unwind(|| fft::<8, _, _>(&input[..7])).is_err());
    assert!(std::panic::catch_unwind(|| ifft::<8, _, _>(&input)).is_err());
}
