#![allow(clippy::excessive_precision)]
#![forbid(unsafe_code)]
// #![no_std]

use core::f64::consts::PI;

use num_complex::Complex;
use trig_const::{cos, cosh, sin, sinh};

const I: Complex<f64> = Complex::new(0.0, 1.0);
const NEG_2: Complex<f64> = Complex::new(-2.0, 0.0);
const TWIDDLE_CONST: Complex<f64> = cmul(cmul(I, NEG_2), Complex::new(PI, 0.0));

// monarch_derive::generate_switch!();
// monarch_derive::generate_powers_of_two!();
// monarch_derive::generate_coprimes!();
// monarch_derive::generate_mixed_radix!();
// monarch_derive::generate_primes!();
// monarch_derive::generate_iffts!();

pub const fn fft<const N: usize>(input: &[Complex<f64>], output: &mut [Complex<f64>])
where
    [(); N / 2]:,
{
    if N == 1 {
        assert!(input.len() > 0);
        assert!(output.len() > 0);

        output[0] = input[0];
    } else if N == 2 {
        assert!(input.len() > 1);
        assert!(output.len() > 1);

        output[0] = cadd(input[0], input[1]);
        output[1] = csub(input[0], input[1]);
    } else if N.is_power_of_two() {
        assert!(input.len() >= N);
        assert!(output.len() >= N);

        let half = N / 2;

        let mut array = [Complex::new(0.0, 0.0); N];
        let mut array_out = [Complex::new(0.0, 0.0); N];
        let (even, odd) = array.split_at_mut(half);
        let (even_out, odd_out) = array_out.split_at_mut(half);
        // Fill evens
        let mut i = 0;
        while i < half {
            even[i] = input[i * 2];
            i += 1;
        }
        // Fill odds
        let mut i = 0;
        while i < half {
            odd[i] = input[i * 2 + 1];
            i += 1;
        }
        fft::<{ N / 2 }>(&even, &mut even_out);

        let mut array_2 = [Complex::new(0.0, 0.0); N];
        let (t, _) = array_2.split_at_mut(half);
        let mut i = 0;
        while i < half {
            let mult_terms = cmul(TWIDDLE_CONST, Complex::new(i as f64, 0.0));
            let pre_exp = cdiv(mult_terms, Complex::new(N as f64, 0.0));
            t[i] = cexp(pre_exp);
            t[i] = cmul(t[i], odd[i]);
            i += 1;
        }

        let mut i = 0;
        while i < half {
            output[i] = cadd(even[i], t[i]);
            i += 1;
        }
        let mut i = 0;
        while i < half {
            output[i] = csub(even[i], t[i]);
            i += 1;
        }
    } else {
        unreachable!()
    }
}

const fn cadd(lhs: Complex<f64>, rhs: Complex<f64>) -> Complex<f64> {
    Complex::new(lhs.re + rhs.re, lhs.im + rhs.im)
}

const fn csub(lhs: Complex<f64>, rhs: Complex<f64>) -> Complex<f64> {
    Complex::new(lhs.re - rhs.re, lhs.im - rhs.im)
}

const fn cmul(lhs: Complex<f64>, rhs: Complex<f64>) -> Complex<f64> {
    Complex::new(
        lhs.re * rhs.re - lhs.im * rhs.im,
        lhs.re * rhs.im + lhs.im * rhs.re,
    )
}

const fn cdiv(lhs: Complex<f64>, rhs: Complex<f64>) -> Complex<f64> {
    let norm_sqr = rhs.re * rhs.re + rhs.im * rhs.im;
    Complex::new(
        (lhs.re * rhs.re + lhs.im * rhs.im) / norm_sqr,
        (lhs.im * rhs.re - lhs.re * rhs.im) / norm_sqr,
    )
}

const fn cexp(x: Complex<f64>) -> Complex<f64> {
    let Complex { re, im } = x;

    if re.is_infinite() {
        if re < 0.0 {
            if !im.is_finite() {
                return Complex::ZERO;
            }
        } else if im == 0.0 || !im.is_finite() {
            if im.is_infinite() {
                return Complex::new(re, f64::NAN);
            }
            return Complex::new(re, im);
        }
    } else if re.is_nan() && im == 0.0 {
        return x;
    }

    cpol(exp(re), im)
}

const fn cpol(r: f64, theta: f64) -> Complex<f64> {
    Complex::new(r * cos(theta), r * sin(theta))
}

const fn csin(x: Complex<f64>) -> Complex<f64> {
    Complex::new(sin(x.re) * cosh(x.im), cos(x.re) * sinh(x.im))
}

const fn ccos(x: Complex<f64>) -> Complex<f64> {
    Complex::new(cos(x.re) * cosh(x.im), -sin(x.re) * sinh(x.im))
}

/// e^x
const fn exp(x: f64) -> f64 {
    let mut i = 1;
    let mut s = 1.0;

    while i < 16 {
        s += expi(x, i) / factorial(i as f64);
        i += 1;
    }

    s
}

/// x^pow
const fn expi(x: f64, mut pow: usize) -> f64 {
    let mut o = 1.0;

    while pow > 0 {
        o *= x;
        pow -= 1;
    }

    o
}

/// Factorial (x!)
const fn factorial(mut x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else {
        let mut s = 1.0;
        while x > 1.0 {
            s *= x;
            x -= 1.0;
        }
        s
    }
}

const fn compute_twiddle(index: usize, fft_len: usize) -> Complex<f64> {
    let constant = -2.0 * PI / fft_len as f64;
    let angle = constant * index as f64;

    Complex::new(cos(angle), sin(angle))
}

#[cfg(test)]
mod rustfft_tests {
    use num_complex::Complex;

    use crate::{cexp, cmul};

    macro_rules! float_eq {
        ($lhs:expr, $rhs:expr) => {
            assert!(($lhs - $rhs).abs() < 0.0001, "lhs: {}, rhs: {}", $lhs, $rhs);
        };
    }

    #[test]
    fn t() {
        let mut plan = rustfft::FftPlannerScalar::<f64>::new();
        let fft = plan.plan_fft_forward(2);
        let mut d = [Complex::new(0.0, 0.0); 2];
        let mut out = d.clone();
        fft.process_with_scratch(&mut d, &mut out);
    }

    #[test]
    fn test_exp() {
        for re in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] {
            for im in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] {
                let c = Complex::new(re, im);

                let stdlib = c.exp();
                let monarch = cexp(c);
                float_eq!(stdlib.re, monarch.re);
                float_eq!(stdlib.im, monarch.im);
            }
        }
    }

    #[test]
    fn test_cmul() {
        for re1 in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] {
            for im1 in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] {
                for re2 in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] {
                    for im2 in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] {
                        let v1 = Complex::new(re1, im1);
                        let v2 = Complex::new(re2, im2);

                        assert_eq!(v1 * v2, cmul(v1, v2));
                    }
                }
            }
        }
    }
}
