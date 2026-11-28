use core::{any::TypeId, ptr};

use num_complex::Complex;
use num_traits::Float;

#[cfg_attr(target_arch = "aarch64", path = "simd/neon.rs")]
#[cfg_attr(target_arch = "x86_64", path = "simd/sse2.rs")]
mod arch;

#[cfg(test)]
mod tests;

#[inline(always)]
pub(super) fn try_fft<const N: usize, const INVERSE: bool, T: Float>(
    input: &[Complex<T>; N],
) -> Option<[Complex<T>; N]> {
    let id = typeid::of::<T>();
    let single = id == TypeId::of::<f32>();
    if !matches!(N, 5 | 8) || (!single && (N != 8 || id != TypeId::of::<f64>())) {
        return None;
    }

    let mut output = [Complex::new(T::zero(), T::zero()); N];
    let src = ptr::from_ref(input);
    let dst = ptr::from_mut(&mut output);
    // SAFETY: The type check identifies the exact primitive (neither
    // contains lifetimes). Each repr(C) complex array has 2 * N contiguous
    // scalars of that type. Input and output are disjoint, and N is 5 or 8.
    unsafe {
        if N == 5 {
            write_fft5::<INVERSE>(src.cast(), dst.cast());
        } else if single {
            write_fft8::<INVERSE, arch::F32Pair>(src.cast(), dst.cast());
        } else {
            write_fft8::<INVERSE, F64Pair>(src.cast(), dst.cast());
        }
    }
    Some(output)
}

trait Pair: Copy {
    type Scalar;

    // Each pointer must cover four scalars. Vector alignment is optional.
    unsafe fn load(src: *const Self::Scalar) -> Self;
    unsafe fn store(self, dst: *mut Self::Scalar);
    fn add(self, rhs: Self) -> Self;
    fn sub(self, rhs: Self) -> Self;
    fn turn<const INVERSE: bool>(self) -> Self;
    fn diagonals<const INVERSE: bool>(self) -> Self;
    fn join<const LEFT: usize, const RIGHT: usize>(self, rhs: Self) -> Self;
}

/// Both pointers must cover ten scalars in disjoint buffers.
#[inline(always)]
unsafe fn write_fft5<const I: bool>(src: *const f32, dst: *mut f32) {
    use arch::F32Pair;

    // SAFETY: All three loads fit the ten input scalars.
    let (base, bc, de) = unsafe {
        (
            F32Pair::load(src),
            F32Pair::load(src.add(2)),
            F32Pair::load(src.add(6)),
        )
    };
    let a = base.join::<0, 0>(base);
    let ed = de.join::<1, 0>(de);
    let sum = bc.add(ed);
    let diff = bc.sub(ed);
    let (p, q) = (sum.join::<0, 0>(sum), sum.join::<1, 1>(sum));
    let (r, s) = (diff.join::<0, 0>(diff), diff.join::<1, 1>(diff));
    let [cos1, cos2, sin1, sin2] = crate::fixed::radix5_coefficients::<f32>();

    // Calculate bins 1 and 2 together, then their conjugate-twiddle partners.
    let even = a.add(p.scale(cos1, cos2)).add(q.scale(cos2, cos1));
    let odd = r.scale(sin1, sin2).add(s.scale(sin2, -sin1)).turn::<I>();
    let upper = even.sub(odd);
    // SAFETY: These stores cover bins 0, 1–2, and 3–4 respectively.
    unsafe {
        a.add(p).add(q).store_first(dst);
        even.add(odd).store(dst.add(2));
        upper.join::<1, 0>(upper).store(dst.add(6));
    }
}

/// Both pointers must cover sixteen scalars in disjoint buffers.
#[inline(always)]
unsafe fn write_fft8<const I: bool, P: Pair>(src: *const P::Scalar, dst: *mut P::Scalar) {
    // SAFETY: Four groups of four scalars fit each caller-provided buffer.
    let [x01, x23, x45, x67] = unsafe {
        [
            P::load(src),
            P::load(src.add(4)),
            P::load(src.add(8)),
            P::load(src.add(12)),
        ]
    };
    let (s01, s23) = (x01.add(x45), x23.add(x67));
    let (d01, d23) = (x01.sub(x45), x23.sub(x67));
    let diagonals = d01.join::<1, 1>(d23).diagonals::<I>();

    // Pair corresponding terms of the even and odd four-point FFTs.
    let a = s01.join::<0, 0>(d01);
    let b = s01.join::<1, 0>(diagonals);
    let c = s23.join::<0, 0>(d23.turn::<I>());
    let d = s23.join::<1, 1>(diagonals);
    let (sum, diff) = (a.add(c), a.sub(c));
    let (odd, skew) = (b.add(d), b.sub(d).turn::<I>());

    // SAFETY: The same four groups cover the output, in bin order.
    unsafe {
        sum.add(odd).store(dst);
        diff.add(skew).store(dst.add(4));
        sum.sub(odd).store(dst.add(8));
        diff.sub(skew).store(dst.add(12));
    }
}

// With two f64 lanes per register, a pair is held in two registers.
#[derive(Clone, Copy)]
struct F64Pair([arch::F64; 2]);

impl Pair for F64Pair {
    type Scalar = f64;

    #[inline(always)]
    unsafe fn load(src: *const f64) -> Self {
        // SAFETY: The caller supplied four readable scalars.
        unsafe { Self([arch::load64(src), arch::load64(src.add(2))]) }
    }

    #[inline(always)]
    unsafe fn store(self, dst: *mut f64) {
        // SAFETY: The caller supplied four writable scalars.
        unsafe {
            arch::store64(dst, self.0[0]);
            arch::store64(dst.add(2), self.0[1]);
        }
    }

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        // SAFETY: The selected target supports these register operations.
        Self(core::array::from_fn(|i| unsafe {
            arch::add64(self.0[i], rhs.0[i])
        }))
    }

    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        // SAFETY: The selected target supports these register operations.
        Self(core::array::from_fn(|i| unsafe {
            arch::sub64(self.0[i], rhs.0[i])
        }))
    }

    #[inline(always)]
    fn turn<const I: bool>(self) -> Self {
        Self(self.0.map(arch::turn64::<I>))
    }

    #[inline(always)]
    fn diagonals<const I: bool>(self) -> Self {
        let [a, b] = self.0.map(arch::scale64);
        // SAFETY: The selected target supports these register operations.
        unsafe {
            Self([
                arch::add64(a, arch::turn64::<I>(a)),
                arch::sub64(arch::turn64::<I>(b), b),
            ])
        }
    }

    #[inline(always)]
    fn join<const L: usize, const R: usize>(self, rhs: Self) -> Self {
        Self([self.0[L], rhs.0[R]])
    }
}
