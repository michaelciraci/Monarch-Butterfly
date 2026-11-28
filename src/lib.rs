#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![no_std]

use num_complex::Complex;
use num_traits::{Float, FloatConst};

mod fixed;

const MAX_FFT_LEN: usize = 200;

/// Computes the forward FFT of exactly `N` samples, for `1 <= N <= 200`.
pub fn fft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    transform::<N, false, _>(input.as_ref())
}

/// Computes the unnormalized inverse FFT of exactly `N` samples, for `1 <= N <= 200`.
pub fn ifft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(
    input: A,
) -> [Complex<T>; N] {
    transform::<N, true, _>(input.as_ref())
}

fn transform<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>],
) -> [Complex<T>; N] {
    assert_eq!(input.len(), N);
    assert!(
        (1..=MAX_FFT_LEN).contains(&N),
        "FFT size must be between 1 and 200"
    );

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
        return radix2_fft::<N, INVERSE, _>(input);
    }

    let zero = Complex::new(T::zero(), T::zero());
    let mut output = [zero; N];
    let mut scratch = [zero; N];
    let roots = core::array::from_fn(|i| twiddle::<INVERSE, T>(Roots::<N>::FULL[i]));
    fft_stage::<N, INVERSE, 0, _>(input, &mut output, &mut scratch, &roots);
    output
}

fn radix2_fft<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>],
) -> [Complex<T>; N] {
    if N == 1 {
        return core::array::from_fn(|i| input[i]);
    }
    let shift = usize::BITS - N.trailing_zeros();
    let mut output = core::array::from_fn(|i| input[i.reverse_bits() >> shift]);

    // In bit-reversed order, adjacent halves form the next butterfly.
    // The supported sizes need at most seven stages. Shorter FFTs skip the rest.
    radix2_stage::<1, INVERSE, _>(&mut output);
    radix2_stage::<2, INVERSE, _>(&mut output);
    radix2_stage::<4, INVERSE, _>(&mut output);
    radix2_stage::<8, INVERSE, _>(&mut output);
    radix2_stage::<16, INVERSE, _>(&mut output);
    radix2_stage::<32, INVERSE, _>(&mut output);
    radix2_stage::<64, INVERSE, _>(&mut output);
    output
}

#[inline(always)]
fn radix2_stage<const HALF: usize, const INVERSE: bool, T: Float>(values: &mut [Complex<T>]) {
    if values.len() < 2 * HALF {
        return;
    }
    let roots: [_; HALF] =
        core::array::from_fn(|k| twiddle::<INVERSE, T>(Roots::<HALF>::HALF_TURN[k]));
    for block in values.chunks_exact_mut(2 * HALF) {
        let (low, high) = block.split_at_mut(HALF);
        for k in 0..HALF {
            let even = low[k];
            let odd = if k == 0 { high[k] } else { high[k] * roots[k] };
            low[k] = even + odd;
            high[k] = even - odd;
        }
    }
}

struct Roots<const N: usize>;

impl<const N: usize> Roots<N> {
    const FULL: [Complex<f64>; N] = make_roots(N);
    const HALF_TURN: [Complex<f64>; N] = make_roots(2 * N);
}

const fn make_roots<const COUNT: usize>(period: usize) -> [Complex<f64>; COUNT] {
    let mut roots = [Complex::new(0.0, 0.0); COUNT];
    let mut k = 0;
    while k < COUNT && period <= MAX_FFT_LEN {
        let angle = core::f64::consts::TAU * k as f64 / period as f64;
        roots[k] = Complex::new(trig_const::cos(angle), trig_const::sin(angle));
        k += 1;
    }
    roots
}

#[inline]
fn twiddle<const INVERSE: bool, T: Float>(root: Complex<f64>) -> Complex<T> {
    let im = if INVERSE { root.im } else { -root.im };
    Complex::new(T::from(root.re).unwrap(), T::from(im).unwrap())
}

const fn split_radix(len: usize) -> usize {
    // The 6-point butterfly combines a pair of 3-point FFTs without twiddles.
    if len % 6 == 0 {
        return 6;
    }

    let mut radix = 2;
    while radix * radix <= len {
        if len % radix == 0 {
            return radix;
        }
        radix += 1;
    }
    len
}

// At depth D, the plan has already split off its first D radices.
struct Stage<const N: usize, const DEPTH: usize>;

impl<const N: usize, const DEPTH: usize> Stage<N, DEPTH> {
    const SHAPE: (usize, usize, usize) = {
        // transform rejects unsupported sizes before entering any FFT stage.
        let mut len = if N > MAX_FFT_LEN { 1 } else { N };
        let mut stride = 1;
        let mut level = 0;
        loop {
            let radix = split_radix(len);
            if level == DEPTH || len <= 1 {
                break (len, radix, stride);
            }
            len /= radix;
            stride *= radix;
            level += 1;
        }
    };
}

// Child stages exchange output and scratch. Only the input offset varies between rows.
fn fft_stage<const N: usize, const INVERSE: bool, const DEPTH: usize, T: Float>(
    input: &[Complex<T>],
    output: &mut [Complex<T>],
    scratch: &mut [Complex<T>],
    roots: &[Complex<T>; N],
) {
    let (len, radix, stride) = Stage::<N, DEPTH>::SHAPE;
    if len <= 1 {
        output[0] = input[0];
        return;
    }
    if len == 2 {
        output[0] = input[0] + input[stride];
        output[1] = input[0] - input[stride];
        return;
    }
    if len == radix && fixed::has_butterfly(radix) {
        fixed::butterfly::<INVERSE, _>(radix, |i| input[i * stride], |i, bin| output[i] = bin);
        return;
    }
    if len == radix {
        // For odd lengths, roots for r and len-r are conjugates. Each pair
        // contributes to both frequency bins k and len-k.
        let half = len / 2;
        let (even, odd) = scratch.split_at_mut(half);
        even.fill(input[0]);
        odd[..half].fill(Complex::new(T::zero(), T::zero()));
        let mut dc = input[0];
        for r in 1..=half {
            let a = input[r * stride];
            let b = input[(len - r) * stride];
            let sum = a + b;
            let diff = a - b;
            dc = dc + sum;
            for k in 0..half {
                let root = roots[(r * (k + 1) % len) * stride];
                even[k] = even[k] + sum * root.re;
                odd[k] = odd[k] + diff * root.im;
            }
        }
        output[0] = dc;
        for k in 0..half {
            let rotated = Complex::new(-odd[k].im, odd[k].re);
            output[k + 1] = even[k] + rotated;
            output[len - 1 - k] = even[k] - rotated;
        }
        return;
    }

    let width = len / radix;
    for r in 0..radix {
        let row = &mut scratch[r * width..(r + 1) * width];
        let work = &mut output[r * width..(r + 1) * width];
        let samples = &input[r * stride..];
        // Stable const generics need literal successors. N <= 200 has at most
        // seven prime factors, so the deepest stage (depth 6) is always a leaf.
        match DEPTH {
            0 => fft_stage::<N, INVERSE, 1, _>(samples, row, work, roots),
            1 => fft_stage::<N, INVERSE, 2, _>(samples, row, work, roots),
            2 => fft_stage::<N, INVERSE, 3, _>(samples, row, work, roots),
            3 => fft_stage::<N, INVERSE, 4, _>(samples, row, work, roots),
            4 => fft_stage::<N, INVERSE, 5, _>(samples, row, work, roots),
            5 => fft_stage::<N, INVERSE, 6, _>(samples, row, work, roots),
            _ => unreachable!(),
        }
    }

    if radix == 2 {
        let (even, odd) = scratch.split_at(width);
        let (low, high) = output.split_at_mut(width);
        for k in 0..width {
            let rotated = roots[k * stride] * odd[k];
            low[k] = even[k] + rotated;
            high[k] = even[k] - rotated;
        }
    } else if fixed::has_butterfly(radix) {
        // Each column is a radix-point FFT after rotation by the stage twiddles.
        for k in 0..width {
            fixed::butterfly::<INVERSE, _>(
                radix,
                |row| {
                    let value = scratch[row * width + k];
                    if row == 0 || k == 0 {
                        value
                    } else {
                        roots[row * k * stride] * value
                    }
                },
                |bin, value| output[bin * width + k] = value,
            );
        }
    } else {
        for (k, bin) in output[..len].iter_mut().enumerate() {
            let col = k % width;
            let mut sum = scratch[col];
            for r in 1..radix {
                let root = roots[(r * k % len) * stride];
                sum = sum + root * scratch[r * width + col];
            }
            *bin = sum;
        }
    }
}
