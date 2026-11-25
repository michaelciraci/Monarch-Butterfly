// #![doc = include_str!("../README.md")]
#![allow(clippy::excessive_precision)]
#![forbid(unsafe_code)]
#![no_std]

use core::{array::from_fn, fmt::Debug};

use num_complex::Complex;
use num_traits::{ConstZero, Float, FloatConst};

/// Computes the forward FFT for lengths whose only prime factors are 2 and 3.
/// A zero-length input returns an empty array.
/// Trigonometric coefficients are evaluated at compile time. Supported lengths
/// through 200 use inlined, fixed-size kernels. Larger transforms look up their
/// twiddles, combining table entries for lengths above 4096.
///
/// # Panics
///
/// Panics if the input length differs from `N` or `N` has another prime factor.
#[inline(always)]
pub fn fft<const N: usize, T: Float + FloatConst + ConstZero + Debug, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    let input = input.as_ref();
    assert_eq!(input.len(), N);
    if N <= 1 {
        return from_fn(|i| input[i]);
    }

    if N <= MAX_KERNEL_LEN {
        return small_fft(input);
    }

    let binary_digits = N.trailing_zeros();
    let binary_len = 1usize << binary_digits;
    let ternary_len = N / binary_len;
    let mut remaining = ternary_len;
    let mut ternary_digits = 0;
    while remaining % 3 == 0 {
        remaining /= 3;
        ternary_digits += 1;
    }
    assert_eq!(remaining, 1, "unsupported FFT length {}", N);

    // Each radix-2 subtransform reads every ternary_len-th input, in bit-reversed
    // order. Order the subtransforms by reversed base-3 digits for the later stages.
    let mut out = from_fn(|idx| {
        let binary_index = if binary_digits == 0 {
            0
        } else {
            (idx % binary_len).reverse_bits() >> (usize::BITS - binary_digits)
        };
        let mut group = idx / binary_len;
        let mut offset = 0;
        for _ in 0..ternary_digits {
            offset = 3 * offset + group % 3;
            group /= 3;
        }
        input[binary_index * ternary_len + offset]
    });

    let mut width = 1;
    while width < binary_len {
        let len = 2 * width;
        let stride = N / len;
        for k in 0..width {
            let w = StageTwiddles::<N>::get::<T>(k * stride);
            for block in out.chunks_exact_mut(len) {
                let even = block[k];
                let odd = w * block[k + width];
                block[k] = even + odd;
                block[k + width] = even - odd;
            }
        }
        width = len;
    }

    while width < N {
        let len = 3 * width;
        let stride = N / len;
        for k in 0..width {
            let index = k * stride;
            let w = StageTwiddles::<N>::get::<T>(index);
            let w2 = StageTwiddles::<N>::get::<T>(2 * index);
            for block in out.chunks_exact_mut(len) {
                let [a, b, c] =
                    radix3_butterfly(block[k], w * block[k + width], w2 * block[k + 2 * width]);
                block[k] = a;
                block[k + width] = b;
                block[k + 2 * width] = c;
            }
        }
        width = len;
    }

    out
}

const fn forward_twiddle(k: usize, period: usize) -> Complex<f64> {
    if k == 0 {
        return Complex::new(1.0, 0.0);
    }
    if period % 4 == 0 {
        if k == period / 4 {
            return Complex::new(0.0, -1.0);
        }
        if k == period - period / 4 {
            return Complex::new(0.0, 1.0);
        }
    }
    if period % 2 == 0 && k == period / 2 {
        return Complex::new(-1.0, 0.0);
    }
    let angle = -core::f64::consts::TAU * k as f64 / period as f64;
    Complex::new(trig_const::cos(angle), trig_const::sin(angle))
}

const fn twiddle_table<const COUNT: usize>(period: usize) -> [Complex<f64>; COUNT] {
    let mut roots = [Complex::new(1.0, 0.0); COUNT];
    let mut k = 1;
    while k < COUNT {
        roots[k] = forward_twiddle(k, period);
        k += 1;
    }
    roots
}

// Bound const-evaluation work and read-only data per transform size.
const TABLE_LEN: usize = 1 << 12;
struct StageTwiddles<const N: usize>;

impl<const N: usize> StageTwiddles<N> {
    const LOW: &'static [Complex<f64>] = if N <= TABLE_LEN {
        &twiddle_table::<N>(N)
    } else {
        &twiddle_table::<TABLE_LEN>(N)
    };

    // W_N^(a+b) = W_N^a * W_N^b. For large N, one constant per high
    // exponent bit extends the low table to cover the entire transform.
    const HIGH: &'static [Complex<f64>; usize::BITS as usize] = &{
        let mut powers = [Complex::new(1.0, 0.0); usize::BITS as usize];
        let mut bit = TABLE_LEN.trailing_zeros();
        while bit < usize::BITS && (1usize << bit) < N {
            powers[bit as usize] = forward_twiddle(1usize << bit, N);
            bit += 1;
        }
        powers
    };

    #[inline(always)]
    fn get<T: Float>(k: usize) -> Complex<T> {
        let mut w = Self::LOW[k % TABLE_LEN];
        if N > TABLE_LEN {
            let mut high = k & !(TABLE_LEN - 1);
            while high != 0 {
                w *= Self::HIGH[high.trailing_zeros() as usize];
                high &= high - 1;
            }
        }
        cast_twiddle(w)
    }
}

#[inline(always)]
fn cast_twiddle<T: Float>(w: Complex<f64>) -> Complex<T> {
    Complex::new(T::from(w.re).unwrap(), T::from(w.im).unwrap())
}

#[inline(always)]
fn radix3_butterfly<T: Float + FloatConst + ConstZero>(
    a: Complex<T>,
    b: Complex<T>,
    c: Complex<T>,
) -> [Complex<T>; 3] {
    let half = T::from(0.5).unwrap();
    // sqrt(3) / 2
    let sqrt3_over_2 = T::from(0.86602540378443864676).unwrap();

    let s = b + c;
    let d = b - c;

    let y0 = a + s;

    let common_re = a.re - half * s.re;
    let common_im = a.im - half * s.im;
    let rot_re = sqrt3_over_2 * d.im;
    let rot_im = sqrt3_over_2 * d.re;

    let y1 = Complex::new(common_re + rot_re, common_im - rot_im);
    let y2 = Complex::new(common_re - rot_re, common_im + rot_im);

    [y0, y1, y2]
}

include!(concat!(env!("OUT_DIR"), "/kernels.rs"));

#[inline(always)]
fn radix2_butterfly<T: Float>(a: Complex<T>, b: Complex<T>) -> [Complex<T>; 2] {
    [a + b, a - b]
}

#[inline(always)]
fn radix4_butterfly<T: Float>(
    a: Complex<T>,
    b: Complex<T>,
    c: Complex<T>,
    d: Complex<T>,
) -> [Complex<T>; 4] {
    let even_sum = a + c;
    let even_diff = a - c;
    let odd_sum = b + d;
    let odd_diff = b - d;
    let rotated = Complex::new(odd_diff.im, -odd_diff.re);
    [
        even_sum + odd_sum,
        even_diff + rotated,
        even_sum - odd_sum,
        even_diff - rotated,
    ]
}

#[inline(always)]
fn rotate<T: Float>(v: Complex<T>, re: T, im: T) -> Complex<T> {
    if cfg!(any(
        all(target_arch = "aarch64", target_feature = "neon"),
        target_feature = "fma"
    )) {
        Complex::new(v.re.mul_add(re, -(v.im * im)), v.im.mul_add(re, v.re * im))
    } else {
        v * Complex::new(re, im)
    }
}
