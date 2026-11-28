//! The NEON target feature is checked at compile time for this module.

use core::arch::aarch64::*;

use super::Pair;

const SIGN32: u32 = 1 << 31;
const SIGN64: u64 = 1 << 63;

#[derive(Clone, Copy)]
pub(super) struct F32Pair(float32x4_t);

impl Pair for F32Pair {
    type Scalar = f32;

    #[inline(always)]
    unsafe fn load(src: *const f32) -> Self {
        // SAFETY: The caller supplied four readable scalars.
        unsafe { Self(vld1q_f32(src)) }
    }

    #[inline(always)]
    unsafe fn store(self, dst: *mut f32) {
        // SAFETY: The caller supplied four writable scalars.
        unsafe { vst1q_f32(dst, self.0) }
    }

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        unsafe { Self(vaddq_f32(self.0, rhs.0)) }
    }

    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        unsafe { Self(vsubq_f32(self.0, rhs.0)) }
    }

    #[inline(always)]
    fn turn<const I: bool>(self) -> Self {
        let signs = if I {
            [SIGN32, 0, SIGN32, 0]
        } else {
            [0, SIGN32, 0, SIGN32]
        };
        // SAFETY: The sign mask has four initialized lanes.
        unsafe {
            let swapped = vreinterpretq_u32_f32(vrev64q_f32(self.0));
            Self(vreinterpretq_f32_u32(veorq_u32(
                swapped,
                vld1q_u32(signs.as_ptr()),
            )))
        }
    }

    #[inline(always)]
    fn diagonals<const I: bool>(self) -> Self {
        unsafe {
            let scaled = Self(vmulq_n_f32(self.0, core::f32::consts::FRAC_1_SQRT_2));
            let signs = vld1q_u32([0, 0, SIGN32, SIGN32].as_ptr());
            let signed = vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(scaled.0), signs));
            Self(signed).add(scaled.turn::<I>())
        }
    }

    #[inline(always)]
    fn join<const L: usize, const R: usize>(self, rhs: Self) -> Self {
        unsafe {
            let a = vreinterpretq_f64_f32(self.0);
            let b = vreinterpretq_f64_f32(rhs.0);
            match (L, R) {
                (0, 0) => Self(vreinterpretq_f32_f64(vzip1q_f64(a, b))),
                (1, 1) => Self(vreinterpretq_f32_f64(vzip2q_f64(a, b))),
                (1, 0) => Self(vcombine_f32(vget_high_f32(self.0), vget_low_f32(rhs.0))),
                (0, 1) => Self(vcombine_f32(vget_low_f32(self.0), vget_high_f32(rhs.0))),
                _ => unreachable!(),
            }
        }
    }
}

pub(super) use core::arch::aarch64::{
    float64x2_t as F64, vaddq_f64 as add64, vld1q_f64 as load64, vst1q_f64 as store64,
    vsubq_f64 as sub64,
};

#[inline(always)]
pub(super) fn turn64<const I: bool>(value: F64) -> F64 {
    let signs = if I { [SIGN64, 0] } else { [0, SIGN64] };
    // SAFETY: The sign mask has two initialized lanes.
    unsafe {
        let swapped = vreinterpretq_u64_f64(vextq_f64::<1>(value, value));
        vreinterpretq_f64_u64(veorq_u64(swapped, vld1q_u64(signs.as_ptr())))
    }
}

#[inline(always)]
pub(super) fn scale64(value: F64) -> F64 {
    unsafe { vmulq_n_f64(value, core::f64::consts::FRAC_1_SQRT_2) }
}
