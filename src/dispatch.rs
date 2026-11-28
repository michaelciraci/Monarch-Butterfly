//! Selects the implementation for a complete transform.

use num_complex::Complex;
use num_traits::Float;

#[cfg(any(
    all(
        target_arch = "aarch64",
        target_endian = "little",
        target_feature = "neon"
    ),
    all(target_arch = "x86_64", target_feature = "sse2")
))]
use crate::simd;
use crate::{fixed, mixed, power2, rader, MAX_FFT_LEN};

#[inline(always)]
pub(super) fn transform<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>],
) -> [Complex<T>; N] {
    #[cfg(any(
        all(
            target_arch = "aarch64",
            target_endian = "little",
            target_feature = "neon"
        ),
        all(target_arch = "x86_64", target_feature = "sse2")
    ))]
    if const { matches!(N, 3 | 5 | 8) } {
        if let Some(output) = input
            .try_into()
            .ok()
            .and_then(simd::try_fft::<N, INVERSE, T>)
        {
            return output;
        }
    }

    if const { matches!(N, 3 | 5) } {
        assert_eq!(input.len(), N);
        let mut output = [Complex::new(T::zero(), T::zero()); N];
        fixed::butterfly::<INVERSE, false, _>(N, |i| input[i], |i, bin| output[i] = bin);
        return output;
    }

    transform_general::<N, INVERSE, T>(input)
}

fn transform_general<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>],
) -> [Complex<T>; N] {
    assert_eq!(input.len(), N);
    assert!(
        (1..=MAX_FFT_LEN).contains(&N),
        "FFT size must be between 1 and 200"
    );

    if N == 9 {
        let bins = fixed::fft9::<INVERSE, _>(core::array::from_fn(|i| input[i]));
        return core::array::from_fn(|i| bins[i]);
    }

    if N == 48 {
        let input = input.try_into().expect("validated 48-point input");
        let output = fixed::transform::<48, 6, INVERSE, _>(input);
        return core::array::from_fn(|i| output[i]);
    }

    if N == 96 {
        let input = input.try_into().expect("validated 96-point input");
        let output = fixed::transform::<96, 12, INVERSE, _>(input);
        return core::array::from_fn(|i| output[i]);
    }

    if N.is_power_of_two() {
        let input = input.try_into().expect("validated power-of-two input");
        return power2::transform::<N, INVERSE, _>(input);
    }

    if const { rader::should_use(N) } {
        let input = input.try_into().expect("validated Rader input");
        return rader::transform::<N, INVERSE, _>(input);
    }

    let zero = Complex::new(T::zero(), T::zero());
    let mut output = [zero; N];
    let mut scratch = [zero; N];
    mixed::fft_into::<N, INVERSE, false, _>(input, &mut output, &mut scratch);
    output
}
