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
monarch_derive::generate_coprimes!();
monarch_derive::generate_mixed_radix!();
monarch_derive::generate_primes!();

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
pub fn fft6<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 6] {
    let n = 6;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let left = fft3([x[0], x[2], x[4]]);
    let right = fft3([x[3], x[5], x[1]]);

    let row1 = fft2([left[0], right[0]]);
    let row2 = fft2([left[1], right[1]]);
    let row3 = fft2([left[2], right[2]]);

    [row1[0], row2[1], row3[0], row1[1], row2[0], row3[1]]
}

#[inline]
pub fn fft9<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 9] {
    let n = 9;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let twiddle1: Complex<T> = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle2: Complex<T> = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle4: Complex<T> = Complex::new(
        T::from(-0.93969262078590832).unwrap(),
        T::from(-0.34202014332566888).unwrap(),
    );

    let first = fft3([x[0], x[3], x[6]]);
    let second = fft3([x[1], x[4], x[7]]);
    let third = fft3([x[2], x[5], x[8]]);

    let row0 = fft3([first[0], second[0], third[0]]);
    let row1 = fft3([first[1], second[1] * twiddle1, third[1] * twiddle2]);
    let row2 = fft3([first[2], second[2] * twiddle2, third[2] * twiddle4]);

    [
        row0[0], row1[0], row2[0], row0[1], row1[1], row2[1], row0[2], row1[2], row2[2],
    ]
}

#[inline]
pub fn fft18<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 18] {
    let n = 18;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let twiddle0 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle1 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle2 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle3 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle4 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle5 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle6 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle7 = Complex::new(
        T::from(0.93969262078590842).unwrap(),
        T::from(-0.34202014332566871).unwrap(),
    );
    let twiddle8 = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle9 = Complex::new(T::from(0.5).unwrap(), T::from(-0.8660254037844386).unwrap());
    let twiddle10 = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle11 = Complex::new(
        T::from(-0.1736481776669303).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle12 = Complex::new(T::from(1.0).unwrap(), T::from(0.0).unwrap());
    let twiddle13 = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle14 = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle15 = Complex::new(
        T::from(-0.5).unwrap(),
        T::from(-0.86602540378443881).unwrap(),
    );
    let twiddle16 = Complex::new(
        T::from(-0.93969262078590832).unwrap(),
        T::from(-0.34202014332566888).unwrap(),
    );
    let twiddle17 = Complex::new(
        T::from(-0.93969262078590842).unwrap(),
        T::from(0.34202014332566866).unwrap(),
    );

    let row0 = fft6([x[0], x[3], x[6], x[9], x[12], x[15]]);
    let row1 = fft6([x[1], x[4], x[7], x[10], x[13], x[16]]);
    let row2 = fft6([x[2], x[5], x[8], x[11], x[14], x[17]]);

    let col0 = fft3([row0[0] * twiddle0, row1[0] * twiddle6, row2[0] * twiddle12]);
    let col1 = fft3([row0[1] * twiddle1, row1[1] * twiddle7, row2[1] * twiddle13]);
    let col2 = fft3([row0[2] * twiddle2, row1[2] * twiddle8, row2[2] * twiddle14]);
    let col3 = fft3([row0[3] * twiddle3, row1[3] * twiddle9, row2[3] * twiddle15]);
    let col4 = fft3([row0[4] * twiddle4, row1[4] * twiddle10, row2[4] * twiddle16]);
    let col5 = fft3([row0[5] * twiddle5, row1[5] * twiddle11, row2[5] * twiddle17]);

    [
        col0[0], col1[0], col2[0], col3[0], col4[0], col5[0], col0[1], col1[1], col2[1], col3[1],
        col4[1], col5[1], col0[2], col1[2], col2[2], col3[2], col4[2], col5[2],
    ]
}

#[inline]
pub fn fft27<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 27] {
    let n = 27;
    let x = input.as_ref();
    assert_eq!(n, x.len());

    let row0 = fft9([x[0], x[3], x[6], x[9], x[12], x[15], x[18], x[21], x[24]]);
    let row1 = fft9([x[1], x[4], x[7], x[10], x[13], x[16], x[19], x[22], x[25]]);
    let row2 = fft9([x[2], x[5], x[8], x[11], x[14], x[17], x[20], x[23], x[26]]);

    let twiddle0 = Complex::new(
        T::from(0.97304487057982381).unwrap(),
        T::from(-0.23061587074244017).unwrap(),
    );
    let twiddle1 = Complex::new(
        T::from(0.89363264032341228).unwrap(),
        T::from(-0.44879918020046217).unwrap(),
    );
    let twiddle2 = Complex::new(
        T::from(0.76604444311897801).unwrap(),
        T::from(-0.64278760968653925).unwrap(),
    );
    let twiddle3 = Complex::new(
        T::from(0.59715859170278618).unwrap(),
        T::from(-0.80212319275504373).unwrap(),
    );
    let twiddle4 = Complex::new(
        T::from(0.3960797660391569).unwrap(),
        T::from(-0.918216106880274).unwrap(),
    );
    let twiddle5 = Complex::new(
        T::from(0.17364817766693041).unwrap(),
        T::from(-0.98480775301220802).unwrap(),
    );
    let twiddle6 = Complex::new(
        T::from(-0.058144828910475774).unwrap(),
        T::from(-0.99830815827126817).unwrap(),
    );
    let twiddle7 = Complex::new(
        T::from(-0.28680323271109021).unwrap(),
        T::from(-0.9579895123154889).unwrap(),
    );
    let twiddle8 = Complex::new(
        T::from(-0.68624163786873338).unwrap(),
        T::from(-0.72737364157304896).unwrap(),
    );
    let twiddle9 = Complex::new(
        T::from(-0.93969262078590832).unwrap(),
        T::from(-0.34202014332566888).unwrap(),
    );
    let twiddle10 = Complex::new(
        T::from(-0.99323835774194302).unwrap(),
        T::from(0.11609291412523012).unwrap(),
    );
    let twiddle11 = Complex::new(
        T::from(-0.83548781141293649).unwrap(),
        T::from(0.54950897807080601).unwrap(),
    );

    let col0 = fft3([row0[0], row1[0], row2[0]]);
    let col1 = fft3([row0[1], row1[1] * twiddle0, row2[1] * twiddle1]);
    let col2 = fft3([row0[2], row1[2] * twiddle1, row2[2] * twiddle3]);
    let col3 = fft3([row0[3], row1[3] * twiddle2, row2[3] * twiddle5]);
    let col4 = fft3([row0[4], row1[4] * twiddle3, row2[4] * twiddle7]);
    let col5 = fft3([row0[5], row1[5] * twiddle4, row2[5] * twiddle8]);
    let col6 = fft3([row0[6], row1[6] * twiddle5, row2[6] * twiddle9]);
    let col7 = fft3([row0[7], row1[7] * twiddle6, row2[7] * twiddle10]);
    let col8 = fft3([row0[8], row1[8] * twiddle7, row2[8] * twiddle11]);

    [
        col0[0], col1[0], col2[0], col3[0], col4[0], col5[0], col6[0], col7[0], col8[0], col0[1],
        col1[1], col2[1], col3[1], col4[1], col5[1], col6[1], col7[1], col8[1], col0[2], col1[2],
        col2[2], col3[2], col4[2], col5[2], col6[2], col7[2], col8[2],
    ]
}


#[cfg(test)]
mod tests {
    use num_complex::Complex;

    use crate::*;

    macro_rules! assert_slice_equal {
        ($lhs:expr,$rhs:expr) => {
            assert_slice_equal!($lhs, $rhs, 0.0001);
        };
        ($lhs:expr,$rhs:expr,$tol:expr) => {
            assert_eq!($lhs.len(), $rhs.len());
            for idx in 0..($lhs.len()) {
                assert!(
                    ($lhs[idx].re - $rhs[idx].re).abs() < $tol,
                    "index {} does not match: ({} - {}).abs() < {}",
                    idx,
                    $lhs[idx].re,
                    $rhs[idx].re,
                    $tol
                );
                assert!(
                    ($lhs[idx].im - $rhs[idx].im).abs() < $tol,
                    "index {} does not match: ({} - {}).abs() < {}",
                    idx,
                    $lhs[idx].im,
                    $rhs[idx].im,
                    $tol
                );
            }
        };
    }

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

        assert_slice_equal!(monarch, buf);
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
    fn test_fft5() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(5);
        let mut buf: Vec<_> = (0..5).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

        let monarch = fft5(&buf);

        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft6() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(6);
        let mut buf: Vec<_> = (0..6).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

        let monarch = fft6(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft7() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(7);
        let mut buf: Vec<_> = (0..7).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

        let monarch = fft7(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_8() {
        let a = butterfly((1..9).map(|i| Complex::new(i as f32, 0.0)).collect());
        assert_eq!(a[0], Complex::new(36.0, 0.0));
        assert_eq!(a[6], Complex::new(-4.0, -4.0));
    }

    #[test]
    fn test_fft9() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(9);
        let mut buf: Vec<_> = (0..9).map(|i| Complex::<f64>::new(i as f64, 0.0)).collect();

        let monarch = fft9(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft10() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(10);
        let mut buf: Vec<_> = (0..10)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft10(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft11() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(11);
        let mut buf: Vec<_> = (0..11)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft11(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft12() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(12);
        let mut buf: Vec<_> = (0..12)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft12(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft13() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(13);
        let mut buf: Vec<_> = (0..13)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft13(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft14() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(14);
        let mut buf: Vec<_> = (0..14)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft14(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft15() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(15);
        let mut buf: Vec<_> = (0..15)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft15(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft17() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(17);
        let mut buf: Vec<_> = (0..17)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft17(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft18() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(18);
        let mut buf: Vec<_> = (0..18)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft18(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft19() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(19);
        let mut buf: Vec<_> = (0..19)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft19(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft20() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(20);
        let mut buf: Vec<_> = (0..20)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft20(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft21() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(21);
        let mut buf: Vec<_> = (0..21)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft21(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft22() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(22);
        let mut buf: Vec<_> = (0..22)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft22(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft23() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(23);
        let mut buf: Vec<_> = (0..23)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft23(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft24() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(24);
        let mut buf: Vec<_> = (0..24)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft24(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft25() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(25);
        let mut buf: Vec<_> = (0..25)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft25(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft26() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(26);
        let mut buf: Vec<_> = (0..26)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft26(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft27() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(27);
        let mut buf: Vec<_> = (0..27)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft27(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft28() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(28);
        let mut buf: Vec<_> = (0..28)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft28(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft29() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(29);
        let mut buf: Vec<_> = (0..29)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft29(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
    }

    #[test]
    fn test_fft30() {
        let mut p = rustfft::FftPlanner::new();
        let plan = p.plan_fft_forward(30);
        let mut buf: Vec<_> = (0..30)
            .map(|i| Complex::<f64>::new(i as f64, 0.0))
            .collect();

        let monarch = fft30(&buf);
        plan.process(&mut buf);

        assert_slice_equal!(monarch, buf);
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

        assert_slice_equal!(ba, bv, 0.001);
    }
}
