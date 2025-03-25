//! Experimental FFT library with an emphasis on runtime performance.
//! This currently only works on powers of two. The FFTs are autogenerated
//! with a procedural macro which hand unrolls all the loops with inlined
//! functions to allow the compiler to maximize the SIMD throughput available
//! on the given CPU.
//!
//! This library will use all SIMD features your CPU has available including AVX512,
//! assuming you compile with those features (`RUSTFLAGS="-C target-cpu=native" cargo build`).
//!
//! This library implements FFTs for both `f32` and `f64` for the following sizes:
//!
//! ```no_compile
//! 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048
//! ```
//!
//! If a larger FFT size is needed, just clone the repo and add the needed
//! sizes to the top of `crates\monarch-derive\src\lib.rs` and larger FFTs
//! will be generated. However, this comes at the cost of a longer compile time.
//!
//! ```
//! use monarch_butterfly::*;
//! use num_complex::Complex;
//!
//! let input: Vec<_> = (0..8).map(|i| Complex::new(i as f32, 0.0)).collect();
//! let output_slice = fft8(&input);
//! let output_vec = fft8(input);
//! ```

#![allow(clippy::excessive_precision)]

use num_complex::Complex;
use num_traits::{Float, FloatConst};

const SQRT_3: f64 = 1.7320508075688772;
const SQRT_3_DIV_2: f64 = SQRT_3 / 2.0;

monarch_derive::generate_powers_of_two!();

fn _compute_twiddle<T: Float + FloatConst>(index: usize, fft_len: usize) -> Complex<T> {
    let constant = T::from(-2.0).unwrap() * T::PI() / T::from(fft_len).unwrap();
    // index * -2PI / fft_len
    let angle = constant * T::from(index).unwrap();

    Complex::new(angle.cos(), angle.sin())
}

#[inline]
pub fn fft3<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 3] {
    let n = 3;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let twiddle: Complex<T> = Complex::new(T::from(-0.5).unwrap(), -T::from(SQRT_3_DIV_2).unwrap());

    let xp = x[1] + x[2];
    let xn = x[1] - x[2];
    let sum = x[0] + xp;

    let temp_a = x[0]
        + Complex {
            re: twiddle.re * xp.re,
            im: twiddle.re * xp.im,
        };
    let temp_b = Complex {
        re: -twiddle.im * xn.im,
        im: twiddle.im * xn.re,
    };

    [sum, temp_a + temp_b, temp_a - temp_b]
}

#[inline]
pub fn fft5<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 5] {
    let n = 5;
    let x = input.as_ref();
    assert_eq!(n, x.len());
    let twiddle1: Complex<T> = Complex::new(
        T::from(0.30901699437494745).unwrap(),
        T::from(-0.95105651629515353).unwrap(),
    );
    let twiddle2: Complex<T> = Complex::new(
        T::from(-0.80901699437494734).unwrap(),
        T::from(-0.58778525229247325).unwrap(),
    );

    let x14p = x[1] + x[4];
    let x14n = x[1] - x[4];
    let x23p = x[2] + x[3];
    let x23n = x[2] - x[3];
    let sum = x[0] + x14p + x23p;
    let b14re_a = x[0].re + twiddle1.re * x14p.re + twiddle2.re * x23p.re;
    let b14re_b = twiddle1.im * x14n.im + twiddle2.im * x23n.im;
    let b23re_a = x[0].re + twiddle2.re * x14p.re + twiddle1.re * x23p.re;
    let b23re_b = twiddle2.im * x14n.im + -twiddle1.im * x23n.im;

    let b14im_a = x[0].im + twiddle1.re * x14p.im + twiddle2.re * x23p.im;
    let b14im_b = twiddle1.im * x14n.re + twiddle2.im * x23n.re;
    let b23im_a = x[0].im + twiddle2.re * x14p.im + twiddle1.re * x23p.im;
    let b23im_b = twiddle2.im * x14n.re + -twiddle1.im * x23n.re;

    let out1re = b14re_a - b14re_b;
    let out1im = b14im_a + b14im_b;
    let out2re = b23re_a - b23re_b;
    let out2im = b23im_a + b23im_b;
    let out3re = b23re_a + b23re_b;
    let out3im = b23im_a - b23im_b;
    let out4re = b14re_a + b14re_b;
    let out4im = b14im_a - b14im_b;

    [
        sum,
        Complex::new(out1re, out1im),
        Complex::new(out2re, out2im),
        Complex::new(out3re, out3im),
        Complex::new(out4re, out4im),
    ]
}

#[cfg(test)]
mod tests {
    use num_complex::Complex;

    use crate::*;

    fn butterfly(x: Vec<Complex<f32>>) -> Vec<Complex<f32>> {
        let n = x.len();
        if n <= 1 {
            x
        } else {
            let even: Vec<_> = butterfly(x.iter().step_by(2).cloned().collect());
            let odd: Vec<_> = butterfly(x.iter().skip(1).step_by(2).cloned().collect());

            let t: Vec<_> = (0..(n / 2))
                .map(|k| {
                    Complex::exp(
                        -2.0 * Complex::<f32>::i() * std::f32::consts::PI * k as f32 / n as f32,
                    ) * odd[k]
                })
                .collect();

            [
                (0..(n / 2)).map(|k| even[k] + t[k]).collect::<Vec<_>>(),
                (0..(n / 2)).map(|k| even[k] - t[k]).collect::<Vec<_>>(),
            ]
            .concat()
        }
    }
    #[test]
    fn test_butterfly_2() {
        assert_eq!(
            fft2([Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)]).to_vec(),
            butterfly(vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)])
        );
    }

    #[test]
    fn test_fft3() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(3);
        let mut buf = vec![
            Complex::<f64>::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
        ];

        let monarch = fft3(&buf);

        plan.process(&mut buf);

        assert_eq!(monarch[0], buf[0]);
        assert!((monarch[1].re - buf[1].re).abs() < 0.0000001);
        assert!((monarch[1].im - buf[1].im).abs() < 0.0000001);
        assert!((monarch[2].re - buf[2].re).abs() < 0.0000001);
        assert!((monarch[2].im - buf[2].im).abs() < 0.0000001);
    }

    #[test]
    fn test_fft5() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(5);
        let mut buf = vec![
            Complex::<f64>::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0),
        ];

        let monarch = fft5(&buf);

        plan.process(&mut buf);

        dbg!(&monarch);
        dbg!(&buf);
        assert_eq!(monarch[0], buf[0]);
        assert!((monarch[1].re - buf[1].re).abs() < 0.0000001);
        assert!((monarch[1].im - buf[1].im).abs() < 0.0000001);
        assert!((monarch[2].re - buf[2].re).abs() < 0.0000001);
        assert!((monarch[2].im - buf[2].im).abs() < 0.0000001);
    }

    #[test]
    fn test_8() {
        let a = butterfly((1..9).map(|i| Complex::new(i as f32, 0.0)).collect());
        assert_eq!(a[0], Complex::new(36.0, 0.0));
        assert_eq!(a[6], Complex::new(-4.0, -4.0));
    }

    #[test]
    fn test_butterfly_4() {
        assert_eq!(
            fft4([
                Complex::new(1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(3.0, 0.0),
                Complex::new(4.0, 0.0)
            ])
            .to_vec(),
            butterfly(vec![
                Complex::new(1.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(3.0, 0.0),
                Complex::new(4.0, 0.0),
            ])
        );
    }

    #[test]
    fn test_butterfly_1024() {
        let v: Vec<_> = (0..1024)
            .map(|i: i32| Complex::new(i as f32, i as f32))
            .collect();
        let a: [Complex<f32>; 1024] = (0..1024)
            .map(|i| Complex::new(i as f32, i as f32))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();

        let ba = fft1024(a);
        let bv = butterfly(v);

        assert_eq!(&ba, &*bv);
    }
}
