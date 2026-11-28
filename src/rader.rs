use core::array::from_fn;

use num_complex::Complex;
use num_traits::Float;

use crate::{mixed::fft_into, roots::Roots, MAX_FFT_LEN};

pub(super) const fn should_use(len: usize) -> bool {
    // For a convolution of length 2 * q with q prime, the mixed-radix path
    // still spends most of its time in the quadratic prime butterfly.
    len >= 29
        && len <= MAX_FFT_LEN
        && len <= u8::MAX as usize + 1
        && is_prime(len)
        && !is_prime((len - 1) / 2)
}

const fn is_prime(len: usize) -> bool {
    if len == 2 {
        return true;
    }
    if len < 2 || len % 2 == 0 {
        return false;
    }
    let mut factor = 3;
    while factor * factor <= len {
        if len % factor == 0 {
            return false;
        }
        factor += 2;
    }
    true
}

struct Rader<const N: usize>;

impl<const N: usize> Rader<N> {
    const LEN: usize = N.saturating_sub(1);

    // List the nonzero indices in the order 1, g, g², ..., for a generator modulo N.
    const PERMUTATION: [u8; N] = {
        let mut indices = [0; N];
        if should_use(N) {
            let mut generator = 2;
            loop {
                assert!(generator < N);
                indices[0] = 1;
                let mut power = 1;
                while power < Self::LEN {
                    indices[power] = (indices[power - 1] as usize * generator % N) as u8;
                    if indices[power] == 1 {
                        break;
                    }
                    power += 1;
                }
                if power == Self::LEN {
                    break;
                }
                generator += 1;
            }
        }
        indices
    };

    // Convolution bin j contributes to FFT bin g^(-j).
    const OUTPUT_BINS: [u8; N] = {
        let mut bins = [0; N];
        if should_use(N) {
            let mut j = 0;
            while j < Self::LEN {
                let bin = Self::PERMUTATION[(Self::LEN - j) % Self::LEN] as usize;
                bins[bin] = j as u8;
                j += 1;
            }
        }
        bins
    };

    // FFT of W_N^(g^(-j)). The scale compensates for the unnormalized inverse
    // used to recover the cyclic convolution.
    const KERNEL: [Complex<f64>; N] = {
        let mut spectrum = [Complex::new(0.0, 0.0); N];
        if should_use(N) {
            let roots = &Roots::<N>::FULL;
            let convolution_roots = &Roots::<N>::CONVOLUTION;
            let mut kernel = spectrum;
            let mut j = 0;
            while j < Self::LEN {
                let index = Self::PERMUTATION[(Self::LEN - j) % Self::LEN] as usize;
                kernel[j] = roots[index];
                j += 1;
            }

            let mut k = 0;
            while k < Self::LEN {
                let mut re = 0.0;
                let mut im = 0.0;
                let mut j = 0;
                while j < Self::LEN {
                    let a = kernel[j];
                    let b = convolution_roots[(k * j) % Self::LEN];
                    // Roots stores positive angles. Both factors here use negative angles.
                    re += a.re * b.re - a.im * b.im;
                    im -= a.re * b.im + a.im * b.re;
                    j += 1;
                }
                spectrum[k] = Complex::new(re / Self::LEN as f64, im / Self::LEN as f64);
                k += 1;
            }
        }
        spectrum
    };
}

#[inline]
pub(super) fn transform<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>; N],
) -> [Complex<T>; N] {
    let zero = Complex::new(T::zero(), T::zero());
    let len = Rader::<N>::LEN;
    let permutation = &Rader::<N>::PERMUTATION;
    let mut samples = [zero; N];
    for (sample, &index) in samples[..len].iter_mut().zip(&permutation[..len]) {
        *sample = input[index as usize];
    }

    let mut spectrum = [zero; N];
    let mut scratch = [zero; N];
    fft_into::<N, false, true, _>(&samples, &mut spectrum, &mut scratch);
    let dc = input[0] + spectrum[0];
    let kernel = &Rader::<N>::KERNEL;
    for (k, bin) in spectrum[..len].iter_mut().enumerate() {
        let factor = if INVERSE {
            // FFT(conj(kernel))[k] is conj(FFT(kernel)[-k]).
            let factor = kernel[if k == 0 { 0 } else { len - k }];
            Complex::new(factor.re, -factor.im)
        } else {
            kernel[k]
        };
        *bin = *bin * Complex::new(T::from(factor.re).unwrap(), T::from(factor.im).unwrap());
    }

    fft_into::<N, true, true, _>(&spectrum, &mut samples, &mut scratch);
    let bins = &Rader::<N>::OUTPUT_BINS;
    from_fn(|k| {
        if k == 0 {
            dc
        } else {
            input[0] + samples[bins[k] as usize]
        }
    })
}
