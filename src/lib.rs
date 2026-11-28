#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![no_std]

use num_complex::Complex;
use num_traits::{Float, FloatConst};

const MAX_LEN: usize = 10;

/// Returns the forward FFT of exactly `N` samples, where `1 <= N <= 10`.
///
/// # Panics
/// Panics if `N` is outside the supported range or the input length is not `N`.
pub fn fft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    transform::<N, false, T>(input.as_ref())
}

/// Returns the unnormalized inverse FFT of exactly `N` samples, where `1 <= N <= 10`.
///
/// # Panics
/// Panics if `N` is outside the supported range or the input length is not `N`.
pub fn ifft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    transform::<N, true, T>(input.as_ref())
}

// Each match arm supplies concrete sizes. `adapt` converts that arm's result
// back to the public return type; validation ensures N is the selected size.
#[inline(always)]
fn transform<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>],
) -> [Complex<T>; N] {
    assert!(
        (1..=MAX_LEN).contains(&N),
        "FFT size must be between 1 and 10"
    );
    assert_eq!(input.len(), N, "input length must match FFT size");

    match N {
        1 => adapt(fft1(core::array::from_fn(|i| input[i]))),
        2 => adapt(fft2(core::array::from_fn(|i| input[i]))),
        3 => adapt(fft3::<INVERSE, _>(core::array::from_fn(|i| input[i]))),
        4 => adapt(fft4::<INVERSE, _>(core::array::from_fn(|i| input[i]))),
        5 => adapt(fft5::<INVERSE, _>(core::array::from_fn(|i| input[i]))),
        6 => adapt(fft6::<INVERSE, _>(core::array::from_fn(|i| input[i]))),
        7 => adapt(prime::<7, INVERSE, _>(core::array::from_fn(|i| input[i]))),
        8 => adapt(fft8::<INVERSE, _>(core::array::from_fn(|i| input[i]))),
        9 => adapt(fft9::<INVERSE, _>(core::array::from_fn(|i| input[i]))),
        10 => adapt(fft10::<INVERSE, _>(core::array::from_fn(|i| input[i]))),
        _ => unreachable!(),
    }
}

#[inline(always)]
fn adapt<const N: usize, const M: usize, T: Float>(bins: [Complex<T>; M]) -> [Complex<T>; N] {
    debug_assert_eq!(N, M);
    core::array::from_fn(|i| bins[i])
}

// Roots are computed at compile time. The direction only changes their sign.
struct Roots<const N: usize>;

impl<const N: usize> Roots<N> {
    const VALUES: [Complex<f64>; N] = {
        let mut roots = [Complex::new(0.0, 0.0); N];
        let mut i = 0;
        while i < N && N <= MAX_LEN {
            let angle = core::f64::consts::TAU * i as f64 / N as f64;
            roots[i] = Complex::new(trig_const::cos(angle), trig_const::sin(angle));
            i += 1;
        }
        roots
    };
}

#[inline(always)]
fn twiddle<const N: usize, const INVERSE: bool, T: Float>(power: usize) -> Complex<T> {
    let root = Roots::<N>::VALUES[power % N];
    let im = if INVERSE { root.im } else { -root.im };
    Complex::new(T::from(root.re).unwrap(), T::from(im).unwrap())
}

#[inline(always)]
fn fft1<T: Float>(input: [Complex<T>; 1]) -> [Complex<T>; 1] {
    input
}

#[inline(always)]
fn fft2<T: Float>([a, b]: [Complex<T>; 2]) -> [Complex<T>; 2] {
    [a + b, a - b]
}

#[inline(always)]
fn fft3<const INVERSE: bool, T: Float>([a, b, c]: [Complex<T>; 3]) -> [Complex<T>; 3] {
    let sum = b + c;
    let difference = b - c;
    let midpoint = a - sum * T::from(0.5).unwrap();
    let rotated = if INVERSE {
        Complex::new(-difference.im, difference.re)
    } else {
        Complex::new(difference.im, -difference.re)
    };
    let skew = rotated * T::from(0.866_025_403_784_438_6).unwrap();
    [a + sum, midpoint + skew, midpoint - skew]
}

// Pair conjugate-frequency inputs so both halves of the spectrum share work.
#[inline(always)]
fn fft5<const INVERSE: bool, T: Float>([a, b, c, d, e]: [Complex<T>; 5]) -> [Complex<T>; 5] {
    let p = b + e;
    let q = c + d;
    let r = b - e;
    let s = c - d;
    let first = Roots::<5>::VALUES[1];
    let second = Roots::<5>::VALUES[2];
    let cos1 = T::from(first.re).unwrap();
    let cos2 = T::from(second.re).unwrap();
    let sin1 = T::from(first.im).unwrap();
    let sin2 = T::from(second.im).unwrap();

    let even1 = a + p * cos1 + q * cos2;
    let even2 = a + p * cos2 + q * cos1;
    let rotate = |x: Complex<T>| {
        if INVERSE {
            Complex::new(-x.im, x.re)
        } else {
            Complex::new(x.im, -x.re)
        }
    };
    let odd1 = rotate(r * sin1 + s * sin2);
    let odd2 = rotate(r * sin2 - s * sin1);
    [
        a + p + q,
        even1 + odd1,
        even2 + odd2,
        even2 - odd2,
        even1 - odd1,
    ]
}

// Pair x[r] with x[N-r]: their roots have the same cosine and opposite
// sine. Each pair contributes to both bins k and N-k. Currently used for 7.
#[inline(always)]
fn prime<const N: usize, const INVERSE: bool, T: Float>(input: [Complex<T>; N]) -> [Complex<T>; N] {
    const { assert!(N >= 3 && N % 2 == 1) };
    let mut output = [Complex::new(T::zero(), T::zero()); N];
    let half = N / 2;
    let mut dc = input[0];
    for r in 1..=half {
        dc = dc + input[r] + input[N - r];
    }
    output[0] = dc;

    for k in 1..=half {
        let mut even = input[0];
        let mut odd = Complex::new(T::zero(), T::zero());
        for r in 1..=half {
            // Positive-angle roots let the same sums serve both directions.
            let root = twiddle::<N, true, T>(r * k);
            even = even + (input[r] + input[N - r]) * root.re;
            odd = odd + (input[r] - input[N - r]) * root.im;
        }
        let skew = if INVERSE {
            Complex::new(-odd.im, odd.re)
        } else {
            Complex::new(odd.im, -odd.re)
        };
        output[k] = even + skew;
        output[N - k] = even - skew;
    }
    output
}

// A and B are independent const parameters because stable Rust cannot use
// N / A as an array length or as a const argument to a child transform.
#[inline(always)]
fn compose<const N: usize, const A: usize, const B: usize, const INVERSE: bool, T: Float>(
    input: [Complex<T>; N],
    fft_a: impl Fn([Complex<T>; A]) -> [Complex<T>; A],
    fft_b: impl Fn([Complex<T>; B]) -> [Complex<T>; B],
) -> [Complex<T>; N] {
    const { assert!(A > 0 && B > 0 && N == A * B) };

    let rows: [[Complex<T>; B]; A] =
        core::array::from_fn(|row| fft_b(core::array::from_fn(|i| input[row + A * i])));
    let columns: [[Complex<T>; A]; B] = core::array::from_fn(|k| {
        fft_a(core::array::from_fn(|row| {
            let value = rows[row][k];
            if row == 0 || k == 0 {
                value
            } else {
                value * twiddle::<N, INVERSE, T>(row * k)
            }
        }))
    });
    core::array::from_fn(|k| columns[k % B][k / B])
}

#[inline(always)]
fn fft4<const I: bool, T: Float>(input: [Complex<T>; 4]) -> [Complex<T>; 4] {
    compose::<4, 2, 2, I, T>(input, fft2, fft2)
}

#[inline(always)]
fn fft6<const I: bool, T: Float>(input: [Complex<T>; 6]) -> [Complex<T>; 6] {
    compose::<6, 2, 3, I, T>(input, fft2, fft3::<I, T>)
}

#[inline(always)]
fn fft8<const I: bool, T: Float>(input: [Complex<T>; 8]) -> [Complex<T>; 8] {
    compose::<8, 2, 4, I, T>(input, fft2, fft4::<I, T>)
}

#[inline(always)]
fn fft9<const I: bool, T: Float>(input: [Complex<T>; 9]) -> [Complex<T>; 9] {
    compose::<9, 3, 3, I, T>(input, fft3::<I, T>, fft3::<I, T>)
}

#[inline(always)]
fn fft10<const I: bool, T: Float>(input: [Complex<T>; 10]) -> [Complex<T>; 10] {
    compose::<10, 2, 5, I, T>(input, fft2, fft5::<I, T>)
}
