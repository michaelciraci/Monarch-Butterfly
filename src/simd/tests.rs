extern crate std;

use core::{array::from_fn, fmt::Debug};
use num_traits::FloatConst;

use super::*;
use crate::{fft, fixed, ifft};

fn scalar_fft<const N: usize, const I: bool, T: Float>(input: [Complex<T>; N]) -> [Complex<T>; N] {
    if N == 3 {
        let bins = fixed::fft3::<I, _>(from_fn(|i| input[i]));
        from_fn(|i| bins[i])
    } else if N == 5 {
        let bins = fixed::fft5::<I, _>(from_fn(|i| input[i]));
        from_fn(|i| bins[i])
    } else {
        let bins = fixed::fft8::<I, _>(from_fn(|i| input[i]));
        from_fn(|i| bins[i])
    }
}

fn check<const N: usize, T: Float + FloatConst + Debug>(input: [Complex<T>; N]) {
    let simd = [
        try_fft::<N, false, _>(&input).unwrap(),
        try_fft::<N, true, _>(&input).unwrap(),
    ];
    let scalar = [
        scalar_fft::<N, false, _>(input),
        scalar_fft::<N, true, _>(input),
    ];
    let public = [fft::<N, _, _>(&input), ifft::<N, _, _>(&input)];
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

fn check_inputs<const N: usize, T: Float + FloatConst + Debug>() {
    let zero = Complex::new(T::zero(), T::zero());
    for index in 0..N {
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
                let mut input = [zero; N];
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
            check::<N, _>(from_fn(|_| {
                Complex::new(sample() * scale, sample() * scale)
            }));
        }
    }
}

#[test]
fn both_precisions_match_scalar_for_impulses_special_values_and_random_inputs() {
    check_inputs::<3, f32>();
    check_inputs::<5, f32>();
    check_inputs::<8, f32>();
    check_inputs::<8, f64>();
}

fn unaligned<const N: usize, T: Float + Debug>(write: unsafe fn(*const T, *mut T)) {
    #[repr(C, align(16))]
    struct Guarded<T, const N: usize>(T, [Complex<T>; N], T);

    let guard = T::one();
    let input: Guarded<T, N> = Guarded(
        guard,
        from_fn(|i| Complex::new(T::from(i).unwrap(), -guard)),
        guard,
    );
    let mut output = Guarded(guard, [Complex::new(T::zero(), T::zero()); N], guard);
    let src = ptr::from_ref(&input.1);
    let dst = ptr::from_mut(&mut output.1);
    assert_ne!(src.addr() % 16, 0);
    assert_ne!(dst.addr() % 16, 0);
    // SAFETY: N is 3, 5, or 8, and each disjoint array has 2 * N scalar lanes.
    unsafe { write(src.cast(), dst.cast()) };
    assert_eq!(output.1, scalar_fft::<N, false, _>(input.1));
    assert_eq!((output.0, output.2), (guard, guard));
}

#[test]
fn unaligned_loads_and_stores_preserve_guards() {
    unaligned::<3, _>(write_fft3::<false>);
    unaligned::<5, _>(write_fft5::<false>);
    unaligned::<8, _>(write_fft8::<false, arch::F32Pair>);
    unaligned::<8, _>(write_fft8::<false, F64Pair>);
}

#[test]
fn dispatch_validates_length() {
    let input = [Complex::new(0.0_f32, 0.0); 9];
    assert!(try_fft::<9, false, _>(&input).is_none());
    assert!(try_fft::<3, false, _>(&[Complex::new(0.0_f64, 0.0); 3]).is_none());
    assert!(try_fft::<5, false, _>(&[Complex::new(0.0_f64, 0.0); 5]).is_none());
    assert!(std::panic::catch_unwind(|| fft::<8, _, _>(&input[..7])).is_err());
    assert!(std::panic::catch_unwind(|| ifft::<8, _, _>(&input)).is_err());
}
