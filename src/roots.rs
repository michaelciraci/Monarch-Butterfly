//! Compile-time roots shared by the transform implementations.

use num_complex::Complex;
use num_traits::Float;

use crate::MAX_FFT_LEN;

pub(super) struct Roots<const N: usize>;

impl<const N: usize> Roots<N> {
    pub(super) const FULL: [Complex<f64>; N] = make_roots(N);
    pub(super) const CONVOLUTION: [Complex<f64>; N] = make_roots(N.saturating_sub(1));
}

const fn make_roots<const COUNT: usize>(period: usize) -> [Complex<f64>; COUNT] {
    let mut roots = [Complex::new(0.0, 0.0); COUNT];
    let mut k = 0;
    while k < COUNT && period != 0 && period <= MAX_FFT_LEN {
        let angle = core::f64::consts::TAU * k as f64 / period as f64;
        roots[k] = Complex::new(trig_const::cos(angle), trig_const::sin(angle));
        k += 1;
    }
    roots
}

#[inline]
pub(super) fn twiddle<const INVERSE: bool, T: Float>(root: Complex<f64>) -> Complex<T> {
    let im = if INVERSE { root.im } else { -root.im };
    Complex::new(T::from(root.re).unwrap(), T::from(im).unwrap())
}
