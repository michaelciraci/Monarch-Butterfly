// #![doc = include_str!("../README.md")]
#![allow(clippy::excessive_precision)]
#![forbid(unsafe_code)]
#![no_std]

use core::{array::from_fn, fmt::Debug};

use num_complex::Complex;
use num_traits::{ConstZero, Float, FloatConst};

#[inline(always)]
pub fn fft<const N: usize, T: Float + FloatConst + ConstZero + Debug, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    let input = input.as_ref();
    assert_eq!(input.len(), N);
    if N == 0 {
        from_fn::<Complex<T>, N, _>(|_| Complex::new(T::zero(), T::zero()))
    } else if N == 1 {
        from_fn(|_| input[0])
    } else if N.is_power_of_two() {
        fft_power_2::<N, T>(&input)
    } else {
        // let even = [input[0]];
        // let odd = [input[1]];
        // let t: [Complex<T>; 1usize] = [Complex::exp(
        //     Complex::<T>::i() * T::from(-2.0).unwrap() * T::PI() * T::from(0usize).unwrap()
        //         / T::from(N).unwrap(),
        // ) * odd[0usize]];
        // let out = [even[0usize] + t[0usize], even[0usize] - t[0usize]];
        // from_fn(|i| out[i])
        // } else {
        todo!()
    }
}

#[inline(always)]
fn fft_power_2<const N: usize, T: Float + FloatConst + ConstZero + Debug>(
    x: &[Complex<T>],
) -> [Complex<T>; N] {
    assert_eq!(N, x.len());
    assert!(N > 0 && N.is_power_of_two());

    let bits = N.trailing_zeros();
    let mut out: [Complex<T>; N] = from_fn(|idx| {
        let rev = idx.reverse_bits() >> (usize::BITS - bits);
        x[rev]
    });

    let neg_two = T::from(-2.0).unwrap();
    let pi = T::PI();
    let i_unit = Complex::<T>::i();

    let mut len = 2;
    while len <= N {
        let half = len / 2;
        let stage_len = T::from(len).unwrap();

        let mut start = 0;
        while start < N {
            let mut k = 0;
            while k < half {
                let angle = neg_two * pi * T::from(k).unwrap() / stage_len;
                let w = (i_unit * angle).exp();

                let even_val = out[start + k];
                let odd_val = w * out[start + k + half];

                out[start + k] = even_val + odd_val;
                out[start + k + half] = even_val - odd_val;

                k += 1;
            }
            start += len;
        }

        len *= 2;
    }

    out
}
