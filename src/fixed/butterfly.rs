//! Adapt fixed-size kernels to gathered samples and scattered bins.

use core::array::from_fn;
use num_complex::Complex;
use num_traits::Float;

use super::kernels::{fft3, fft5, fft6};

pub(crate) const fn has_butterfly(radix: usize) -> bool {
    matches!(radix, 3 | 5 | 6)
}

#[inline(always)]
pub(crate) fn butterfly<const INVERSE: bool, const LEAF: bool, T: Float>(
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
