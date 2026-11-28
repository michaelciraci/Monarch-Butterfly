#![doc = include_str!("../README.md")]
#![deny(unsafe_code, unsafe_op_in_unsafe_fn)]
#![no_std]

use num_complex::Complex;
use num_traits::{Float, FloatConst};

mod dispatch;
mod fixed;
mod mixed;
mod power2;
mod rader;
mod roots;
#[cfg(any(
    all(
        target_arch = "aarch64",
        target_endian = "little",
        target_feature = "neon"
    ),
    all(target_arch = "x86_64", target_feature = "sse2")
))]
#[allow(unsafe_code)]
mod simd;

const MAX_FFT_LEN: usize = 200;

/// Computes the forward FFT of exactly `N` samples, for `1 <= N <= 200`.
pub fn fft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    dispatch::transform::<N, false, _>(input.as_ref())
}

/// Computes the unnormalized inverse FFT of exactly `N` samples, for `1 <= N <= 200`.
pub fn ifft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    dispatch::transform::<N, true, _>(input.as_ref())
}
