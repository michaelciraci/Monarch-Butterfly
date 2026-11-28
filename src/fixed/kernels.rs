//! Reusable scalar butterflies. Inputs and outputs are in natural order.

use num_complex::Complex;
use num_traits::Float;

use crate::roots::Roots;

#[inline(always)]
pub(crate) fn fft6<const INVERSE: bool, T: Float>(x: [Complex<T>; 6]) -> [Complex<T>; 6] {
    // A Good–Thomas split into two 3-point FFTs avoids inter-stage twiddles.
    let a = fft3::<INVERSE, _>([x[0], x[2], x[4]]);
    let b = fft3::<INVERSE, _>([x[3], x[5], x[1]]);
    [
        a[0] + b[0],
        a[1] - b[1],
        a[2] + b[2],
        a[0] - b[0],
        a[1] + b[1],
        a[2] - b[2],
    ]
}

#[inline(always)]
pub(crate) fn fft8<const INVERSE: bool, T: Float>(input: [Complex<T>; 8]) -> [Complex<T>; 8] {
    let [a, b, c, d, e, f, g, h] = input;
    let even = fft4::<INVERSE, _>([a + e, b + f, c + g, d + h]);
    // The eighth roots at bins 1 and 3 have equal-magnitude components.
    let scale = T::from(core::f64::consts::FRAC_1_SQRT_2).unwrap();
    let odd_one = (b - f) * scale;
    let odd_three = (d - h) * scale;
    let odd = fft4::<INVERSE, _>([
        a - e,
        odd_one + quarter_turn::<INVERSE, _>(odd_one),
        quarter_turn::<INVERSE, _>(c - g),
        quarter_turn::<INVERSE, _>(odd_three) - odd_three,
    ]);
    [
        even[0], odd[0], even[1], odd[1], even[2], odd[2], even[3], odd[3],
    ]
}

#[inline(always)]
pub(crate) fn fft12<const INVERSE: bool, T: Float>(x: [Complex<T>; 12]) -> [Complex<T>; 12] {
    // The coprime 3-by-4 split uses Good–Thomas permutations to skip twiddles.
    let a = fft3::<INVERSE, _>([x[0], x[4], x[8]]);
    let b = fft3::<INVERSE, _>([x[3], x[7], x[11]]);
    let c = fft3::<INVERSE, _>([x[6], x[10], x[2]]);
    let d = fft3::<INVERSE, _>([x[9], x[1], x[5]]);
    let u = fft4::<INVERSE, _>([a[0], b[0], c[0], d[0]]);
    let v = fft4::<INVERSE, _>([a[1], b[1], c[1], d[1]]);
    let w = fft4::<INVERSE, _>([a[2], b[2], c[2], d[2]]);
    [
        u[0], v[1], w[2], u[3], v[0], w[1], u[2], v[3], w[0], u[1], v[2], w[3],
    ]
}

#[inline(always)]
pub(crate) fn fft3<const INVERSE: bool, T: Float>(input: [Complex<T>; 3]) -> [Complex<T>; 3] {
    let [a, b, c] = input;
    let sum = b + c;
    let midpoint = a - sum * T::from(0.5).unwrap();
    let skew = quarter_turn::<INVERSE, _>(b - c) * T::from(0.866_025_403_784_438_6).unwrap();
    [a + sum, midpoint + skew, midpoint - skew]
}

/// A 3-by-3 Cooley–Tukey transform in natural order.
#[inline(always)]
pub(crate) fn fft9<const INVERSE: bool, T: Float>(x: [Complex<T>; 9]) -> [Complex<T>; 9] {
    use crate::roots::twiddle;

    let a = fft3::<INVERSE, _>([x[0], x[3], x[6]]);
    let b = fft3::<INVERSE, _>([x[1], x[4], x[7]]);
    let c = fft3::<INVERSE, _>([x[2], x[5], x[8]]);
    let roots = &Roots::<9>::FULL;
    let zero = fft3::<INVERSE, _>([a[0], b[0], c[0]]);
    let one = fft3::<INVERSE, _>([
        a[1],
        b[1] * twiddle::<INVERSE, T>(roots[1]),
        c[1] * twiddle::<INVERSE, T>(roots[2]),
    ]);
    let two = fft3::<INVERSE, _>([
        a[2],
        b[2] * twiddle::<INVERSE, T>(roots[2]),
        c[2] * twiddle::<INVERSE, T>(roots[4]),
    ]);
    [
        zero[0], one[0], two[0], zero[1], one[1], two[1], zero[2], one[2], two[2],
    ]
}

#[inline(always)]
pub(crate) fn fft5<const INVERSE: bool, T: Float>(
    [a, b, c, d, e]: [Complex<T>; 5],
) -> [Complex<T>; 5] {
    let p = b + e;
    let q = c + d;
    let r = b - e;
    let s = c - d;
    let [cos1, cos2, sin1, sin2] = radix5_coefficients::<T>();

    // Conjugate twiddles share their cosine terms in bins k and 5 - k.
    let even1 = a + p * cos1 + q * cos2;
    let even2 = a + p * cos2 + q * cos1;
    let odd1 = quarter_turn::<INVERSE, _>(r * sin1 + s * sin2);
    let odd2 = quarter_turn::<INVERSE, _>(r * sin2 - s * sin1);
    [
        a + p + q,
        even1 + odd1,
        even2 + odd2,
        even2 - odd2,
        even1 - odd1,
    ]
}

#[inline(always)]
pub(crate) fn radix5_coefficients<T: Float>() -> [T; 4] {
    let w1 = Roots::<5>::FULL[1];
    let w2 = Roots::<5>::FULL[2];
    [w1.re, w2.re, w1.im, w2.im].map(|x| T::from(x).unwrap())
}

#[inline(always)]
pub(crate) fn fft4<const INVERSE: bool, T: Float>(input: [Complex<T>; 4]) -> [Complex<T>; 4] {
    let [a, b, c, d] = input;
    let even_sum = a + c;
    let even_diff = a - c;
    let odd_sum = b + d;
    let rotated = quarter_turn::<INVERSE, _>(b - d);
    [
        even_sum + odd_sum,
        even_diff + rotated,
        even_sum - odd_sum,
        even_diff - rotated,
    ]
}

#[inline(always)]
pub(crate) fn quarter_turn<const INVERSE: bool, T: Float>(x: Complex<T>) -> Complex<T> {
    if INVERSE {
        Complex::new(-x.im, x.re)
    } else {
        Complex::new(x.im, -x.re)
    }
}
