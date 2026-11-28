//! Fixed-size transforms built from smaller kernels.

use core::array::from_fn;
use num_complex::Complex;
use num_traits::Float;

use super::kernels::{fft12, fft4, fft6, fft8, quarter_turn};
use crate::roots::{twiddle, Roots};

const COLS: usize = 8;

#[inline(always)]
pub(crate) fn transform<const N: usize, const ROWS: usize, const INVERSE: bool, T: Float>(
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
pub(crate) fn rotate<const N: usize, const INVERSE: bool, T: Float>(
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
