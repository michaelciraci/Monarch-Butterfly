//! The SSE2 target feature is checked at compile time for this module.

use core::arch::x86_64::*;

use super::Pair;

#[derive(Clone, Copy)]
pub(super) struct F32Pair(__m128);

impl F32Pair {
    #[inline(always)]
    pub(super) unsafe fn store_first(self, dst: *mut f32) {
        // SAFETY: The caller supplied two writable scalars (eight bytes).
        unsafe { _mm_storel_epi64(dst.cast(), _mm_castps_si128(self.0)) }
    }

    #[inline(always)]
    pub(super) fn scale(self, left: f32, right: f32) -> Self {
        unsafe { Self(_mm_mul_ps(self.0, _mm_setr_ps(left, left, right, right))) }
    }
}

impl Pair for F32Pair {
    type Scalar = f32;

    #[inline(always)]
    unsafe fn load(src: *const f32) -> Self {
        // SAFETY: The caller supplied four readable scalars.
        unsafe { Self(_mm_loadu_ps(src)) }
    }

    #[inline(always)]
    unsafe fn store(self, dst: *mut f32) {
        // SAFETY: The caller supplied four writable scalars.
        unsafe { _mm_storeu_ps(dst, self.0) }
    }

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        unsafe { Self(_mm_add_ps(self.0, rhs.0)) }
    }

    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        unsafe { Self(_mm_sub_ps(self.0, rhs.0)) }
    }

    #[inline(always)]
    fn turn<const I: bool>(self) -> Self {
        unsafe {
            let signs = if I {
                _mm_setr_ps(-0.0, 0.0, -0.0, 0.0)
            } else {
                _mm_setr_ps(0.0, -0.0, 0.0, -0.0)
            };
            Self(_mm_xor_ps(_mm_shuffle_ps::<0xb1>(self.0, self.0), signs))
        }
    }

    #[inline(always)]
    fn diagonals<const I: bool>(self) -> Self {
        unsafe {
            let scaled = Self(_mm_mul_ps(
                self.0,
                _mm_set1_ps(core::f32::consts::FRAC_1_SQRT_2),
            ));
            let signed = _mm_xor_ps(scaled.0, _mm_setr_ps(0.0, 0.0, -0.0, -0.0));
            Self(signed).add(scaled.turn::<I>())
        }
    }

    #[inline(always)]
    fn join<const L: usize, const R: usize>(self, rhs: Self) -> Self {
        unsafe {
            match (L, R) {
                (0, 0) => Self(_mm_movelh_ps(self.0, rhs.0)),
                (1, 1) => Self(_mm_movehl_ps(rhs.0, self.0)),
                (1, 0) => Self(_mm_shuffle_ps::<0x4e>(self.0, rhs.0)),
                (0, 1) => Self(_mm_shuffle_ps::<0xe4>(self.0, rhs.0)),
                _ => unreachable!(),
            }
        }
    }
}

pub(super) use core::arch::x86_64::{
    __m128d as F64, _mm_add_pd as add64, _mm_loadu_pd as load64, _mm_storeu_pd as store64,
    _mm_sub_pd as sub64,
};

#[inline(always)]
pub(super) fn turn64<const I: bool>(value: F64) -> F64 {
    unsafe {
        let signs = if I {
            _mm_setr_pd(-0.0, 0.0)
        } else {
            _mm_setr_pd(0.0, -0.0)
        };
        _mm_xor_pd(_mm_shuffle_pd::<1>(value, value), signs)
    }
}

#[inline(always)]
pub(super) fn scale64(value: F64) -> F64 {
    unsafe { _mm_mul_pd(value, _mm_set1_pd(core::f64::consts::FRAC_1_SQRT_2)) }
}
