use core::array::from_fn;

use num_complex::Complex;
use num_traits::Float;

use crate::fixed;

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
            let bins = radix4::<64, 16, INVERSE, _>(input, fft16::<INVERSE, T>);
            from_fn(|i| bins[i])
        }
        128 => {
            let input = input.as_slice().try_into().expect("128-point input");
            let bins = fixed::transform::<128, 16, INVERSE, _>(input);
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
        16 => output.copy_from_slice(&fft16::<INVERSE, _>(from_fn(sample))),
        32 => output.copy_from_slice(&fixed::transform::<32, 4, INVERSE, _>(&from_fn(sample))),
        _ => unreachable!(),
    }
}

#[inline(always)]
pub(super) fn fft16<const INVERSE: bool, T: Float>(input: [Complex<T>; 16]) -> [Complex<T>; 16] {
    radix4::<16, 4, INVERSE, _>(&input, fixed::fft4::<INVERSE, T>)
}

#[inline(always)]
fn radix4<const N: usize, const Q: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>; N],
    subtransform: impl Fn([Complex<T>; Q]) -> [Complex<T>; Q],
) -> [Complex<T>; N] {
    // Transform the four strided subsequences, then combine each frequency bin.
    let a = subtransform(from_fn(|k| input[4 * k]));
    let b = subtransform(from_fn(|k| input[4 * k + 1]));
    let c = subtransform(from_fn(|k| input[4 * k + 2]));
    let d = subtransform(from_fn(|k| input[4 * k + 3]));
    let mut output = [Complex::new(T::zero(), T::zero()); N];

    for k in 0..Q {
        let [y0, y1, y2, y3] = fixed::fft4::<INVERSE, _>([
            a[k],
            fixed::rotate::<N, INVERSE, _>(b[k], k),
            fixed::rotate::<N, INVERSE, _>(c[k], 2 * k),
            fixed::rotate::<N, INVERSE, _>(d[k], 3 * k),
        ]);
        output[k] = y0;
        output[k + Q] = y1;
        output[k + 2 * Q] = y2;
        output[k + 3 * Q] = y3;
    }
    output
}
