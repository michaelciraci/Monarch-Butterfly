use core::array::from_fn;

use num_complex::Complex;
use num_traits::Float;

use crate::{twiddle, Roots};

const COLS: usize = 8;

pub(super) const fn has_butterfly(radix: usize) -> bool {
    matches!(radix, 3 | 5 | 6)
}

#[inline(always)]
pub(super) fn butterfly<const INVERSE: bool, const LEAF: bool, T: Float>(
    radix: usize,
    mut sample: impl FnMut(usize) -> Complex<T>,
    store: impl FnMut(usize, Complex<T>),
) {
    fn write<const N: usize, T: Copy>(bins: [T; N], mut store: impl FnMut(usize, T)) {
        for (i, bin) in bins.into_iter().enumerate() {
            store(i, bin);
        }
    }

    match radix {
        3 => write(fft3::<INVERSE, _>(from_fn(&mut sample)), store),
        5 => {
            let samples = from_fn(&mut sample);
            #[cfg(any(
                all(
                    target_arch = "aarch64",
                    target_endian = "little",
                    target_feature = "neon"
                ),
                all(target_arch = "x86_64", target_feature = "sse2")
            ))]
            if LEAF {
                if let Some(bins) = crate::simd::try_fft::<5, INVERSE, _>(&samples) {
                    write(bins, store);
                    return;
                }
            }
            write(fft5::<INVERSE, _>(samples), store);
        }
        6 => write(fft6::<INVERSE, _>(from_fn(&mut sample)), store),
        _ => unreachable!(),
    }
}

#[inline(always)]
pub(super) fn transform<const N: usize, const ROWS: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>; N],
) -> [Complex<T>; N] {
    const { assert!(N == ROWS * COLS) };

    // Transform each strided group of eight, then collect matching bins by column.
    let mut columns = [[Complex::new(T::zero(), T::zero()); ROWS]; COLS];
    for row in 0..ROWS {
        let samples = from_fn(|i| input[row + i * ROWS]);
        let bins = fft8::<INVERSE, _>(samples);
        for (column, bin) in columns.iter_mut().zip(bins) {
            column[row] = bin;
        }
    }

    // A const column index specializes the twiddles for each second-stage FFT.
    let bins = [
        finish_column::<0, N, ROWS, INVERSE, _>(&columns[0]),
        finish_column::<1, N, ROWS, INVERSE, _>(&columns[1]),
        finish_column::<2, N, ROWS, INVERSE, _>(&columns[2]),
        finish_column::<3, N, ROWS, INVERSE, _>(&columns[3]),
        finish_column::<4, N, ROWS, INVERSE, _>(&columns[4]),
        finish_column::<5, N, ROWS, INVERSE, _>(&columns[5]),
        finish_column::<6, N, ROWS, INVERSE, _>(&columns[6]),
        finish_column::<7, N, ROWS, INVERSE, _>(&columns[7]),
    ];
    from_fn(|i| bins[i % COLS][i / COLS])
}

#[inline(always)]
fn finish_column<
    const COLUMN: usize,
    const N: usize,
    const ROWS: usize,
    const INVERSE: bool,
    T: Float,
>(
    column: &[Complex<T>; ROWS],
) -> [Complex<T>; ROWS] {
    // Literal row indices let the compiler fold each rotation before the butterfly.
    let rotate = |row: usize| multiply_twiddle::<N, INVERSE, _>(column[row], COLUMN * row);
    match ROWS {
        4 => {
            let bins = fft4::<INVERSE, _>([column[0], rotate(1), rotate(2), rotate(3)]);
            from_fn(|i| bins[i])
        }
        6 => {
            let bins = fft6::<INVERSE, _>([
                column[0],
                rotate(1),
                rotate(2),
                rotate(3),
                rotate(4),
                rotate(5),
            ]);
            from_fn(|i| bins[i])
        }
        12 => {
            let bins = fft12::<INVERSE, _>([
                column[0],
                rotate(1),
                rotate(2),
                rotate(3),
                rotate(4),
                rotate(5),
                rotate(6),
                rotate(7),
                rotate(8),
                rotate(9),
                rotate(10),
                rotate(11),
            ]);
            from_fn(|i| bins[i])
        }
        _ => unreachable!(),
    }
}

#[inline(always)]
fn multiply_twiddle<const N: usize, const INVERSE: bool, T: Float>(
    value: Complex<T>,
    power: usize,
) -> Complex<T> {
    if power == 0 {
        return value;
    }
    if N % 4 == 0 && power == N / 4 {
        return quarter_turn::<INVERSE, _>(value);
    }
    if N % 2 == 0 && power == N / 2 {
        return -value;
    }

    rotate::<N, INVERSE, T>(value, power)
}

#[inline(always)]
pub(super) fn rotate<const N: usize, const INVERSE: bool, T: Float>(
    value: Complex<T>,
    power: usize,
) -> Complex<T> {
    let root = twiddle::<INVERSE, T>(Roots::<N>::FULL[power]);
    if cfg!(any(
        all(target_arch = "aarch64", target_feature = "neon"),
        target_feature = "fma"
    )) {
        Complex::new(
            value.re.mul_add(root.re, -(value.im * root.im)),
            value.im.mul_add(root.re, value.re * root.im),
        )
    } else {
        value * root
    }
}

#[inline(always)]
fn fft6<const INVERSE: bool, T: Float>(x: [Complex<T>; 6]) -> [Complex<T>; 6] {
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
pub(super) fn fft8<const INVERSE: bool, T: Float>(input: [Complex<T>; 8]) -> [Complex<T>; 8] {
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
fn fft12<const INVERSE: bool, T: Float>(x: [Complex<T>; 12]) -> [Complex<T>; 12] {
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
fn fft3<const INVERSE: bool, T: Float>(input: [Complex<T>; 3]) -> [Complex<T>; 3] {
    let [a, b, c] = input;
    let sum = b + c;
    let midpoint = a - sum * T::from(0.5).unwrap();
    let skew = quarter_turn::<INVERSE, _>(b - c) * T::from(0.866_025_403_784_438_6).unwrap();
    [a + sum, midpoint + skew, midpoint - skew]
}

#[inline(always)]
pub(super) fn fft5<const INVERSE: bool, T: Float>(
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
pub(super) fn radix5_coefficients<T: Float>() -> [T; 4] {
    let w1 = Roots::<5>::FULL[1];
    let w2 = Roots::<5>::FULL[2];
    [w1.re, w2.re, w1.im, w2.im].map(|x| T::from(x).unwrap())
}

#[inline(always)]
pub(super) fn fft4<const INVERSE: bool, T: Float>(input: [Complex<T>; 4]) -> [Complex<T>; 4] {
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
pub(super) fn quarter_turn<const INVERSE: bool, T: Float>(x: Complex<T>) -> Complex<T> {
    if INVERSE {
        Complex::new(-x.im, x.re)
    } else {
        Complex::new(x.im, -x.re)
    }
}
