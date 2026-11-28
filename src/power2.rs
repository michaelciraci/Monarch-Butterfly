use core::array::from_fn;

use num_complex::Complex;
use num_traits::Float;

use crate::{
    fixed,
    roots::{twiddle, Roots},
};

#[inline]
pub(super) fn transform<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>; N],
) -> [Complex<T>; N] {
    match N {
        1 => *input,
        2 => {
            let a = input[0];
            let b = input[1];
            from_fn(|i| if i == 0 { a + b } else { a - b })
        }
        4 | 8 | 16 | 32 => {
            let mut output = [Complex::new(T::zero(), T::zero()); N];
            transform_small::<INVERSE, _>(N, |i| input[i], &mut output);
            output
        }
        64 => {
            let input = input.as_slice().try_into().expect("64-point input");
            let bins = radix4::<64, 16, INVERSE, true, _>(input, fft16::<INVERSE, T>);
            from_fn(|i| bins[i])
        }
        128 => {
            let input = input.as_slice().try_into().expect("128-point input");
            // Keep the four rows contiguous through the twiddle multiply and merge.
            let bins = radix4::<128, 32, INVERSE, false, _>(input, |row| {
                radix4::<32, 8, INVERSE, false, _>(&row, fixed::fft8::<INVERSE, T>)
            });
            from_fn(|i| bins[i])
        }
        _ => unreachable!(),
    }
}

// The mixed-radix planner can pass a strided gather to the same kernels.
#[inline(always)]
pub(super) fn transform_small<const INVERSE: bool, T: Float>(
    len: usize,
    sample: impl FnMut(usize) -> Complex<T>,
    output: &mut [Complex<T>],
) {
    let output = &mut output[..len];
    match len {
        4 => output.copy_from_slice(&fixed::fft4::<INVERSE, _>(from_fn(sample))),
        8 => output.copy_from_slice(&fixed::fft8::<INVERSE, _>(from_fn(sample))),
        16 => output.copy_from_slice(&fft16_exact::<INVERSE, _>(from_fn(sample))),
        32 => output.copy_from_slice(&fixed::transform::<32, 4, INVERSE, _>(&from_fn(sample))),
        _ => unreachable!(),
    }
}

#[inline(always)]
fn fft16_exact<const INVERSE: bool, T: Float>(x: [Complex<T>; 16]) -> [Complex<T>; 16] {
    let a = fixed::fft4::<INVERSE, _>([x[0], x[4], x[8], x[12]]);
    let b = fixed::fft4::<INVERSE, _>([x[1], x[5], x[9], x[13]]);
    let c = fixed::fft4::<INVERSE, _>([x[2], x[6], x[10], x[14]]);
    let d = fixed::fft4::<INVERSE, _>([x[3], x[7], x[11], x[15]]);

    // Quarter-turn and diagonal rotations need fewer operations than complex products.
    let zero = fixed::fft4::<INVERSE, _>([a[0], b[0], c[0], d[0]]);
    let one = fixed::fft4::<INVERSE, _>([
        a[1],
        fixed::rotate::<16, INVERSE, _>(b[1], 1),
        diagonal::<INVERSE, _>(c[1]),
        fixed::rotate::<16, INVERSE, _>(d[1], 3),
    ]);
    let two = fixed::fft4::<INVERSE, _>([
        a[2],
        diagonal::<INVERSE, _>(b[2]),
        fixed::quarter_turn::<INVERSE, _>(c[2]),
        fixed::quarter_turn::<INVERSE, _>(diagonal::<INVERSE, _>(d[2])),
    ]);
    let three = fixed::fft4::<INVERSE, _>([
        a[3],
        fixed::rotate::<16, INVERSE, _>(b[3], 3),
        fixed::quarter_turn::<INVERSE, _>(diagonal::<INVERSE, _>(c[3])),
        -fixed::rotate::<16, INVERSE, _>(d[3], 1),
    ]);
    [
        zero[0], one[0], two[0], three[0], zero[1], one[1], two[1], three[1], zero[2], one[2],
        two[2], three[2], zero[3], one[3], two[3], three[3],
    ]
}

#[inline(always)]
fn diagonal<const INVERSE: bool, T: Float>(value: Complex<T>) -> Complex<T> {
    let scale = T::from(core::f64::consts::FRAC_1_SQRT_2).unwrap();
    let re = value.re * scale;
    let im = value.im * scale;
    if INVERSE {
        Complex::new(re - im, re + im)
    } else {
        Complex::new(re + im, im - re)
    }
}

#[inline(always)]
fn fft16<const INVERSE: bool, T: Float>(input: [Complex<T>; 16]) -> [Complex<T>; 16] {
    radix4::<16, 4, INVERSE, true, _>(&input, fixed::fft4::<INVERSE, T>)
}

#[inline(always)]
fn radix4<const N: usize, const Q: usize, const INVERSE: bool, const FUSED: bool, T: Float>(
    input: &[Complex<T>; N],
    subtransform: impl Fn([Complex<T>; Q]) -> [Complex<T>; Q],
) -> [Complex<T>; N] {
    // Transform the four strided subsequences, then combine each frequency bin.
    let a = subtransform(from_fn(|k| input[4 * k]));
    let b = subtransform(from_fn(|k| input[4 * k + 1]));
    let c = subtransform(from_fn(|k| input[4 * k + 2]));
    let d = subtransform(from_fn(|k| input[4 * k + 3]));
    let mut output = [Complex::new(T::zero(), T::zero()); N];
    // Plain complex products give LLVM more room to vectorize the larger stages.
    let rotate = |value, power: usize| {
        if FUSED {
            fixed::rotate::<N, INVERSE, T>(value, power)
        } else {
            value * twiddle::<INVERSE, T>(Roots::<N>::FULL[power])
        }
    };
    for k in 0..Q {
        let [y0, y1, y2, y3] = fixed::fft4::<INVERSE, _>([
            a[k],
            rotate(b[k], k),
            rotate(c[k], 2 * k),
            rotate(d[k], 3 * k),
        ]);
        output[k] = y0;
        output[k + Q] = y1;
        output[k + 2 * Q] = y2;
        output[k + 3 * Q] = y3;
    }
    output
}
