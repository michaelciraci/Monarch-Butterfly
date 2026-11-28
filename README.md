[![Build](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml/badge.svg)](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml)
[![](https://img.shields.io/crates/v/monarch-butterfly)](https://crates.io/crates/monarch-butterfly)
[![](https://docs.rs/monarch-butterfly/badge.svg)](https://docs.rs/monarch-butterfly/)

# Monarch Butterfly

Monarch Butterfly computes forward and inverse FFTs of 1 through 200 complex
samples in `f32` or `f64`. The transform size is a const generic. Power-of-two
sizes use fixed-array kernels that can also serve as mixed-radix
subtransforms. Mixed-radix stage lengths, radices, and strides are derived
from `N` at compile time. The FFT twiddle factors are also computed at compile
time. Selected prime lengths use Rader's algorithm with compile-time
permutations and convolution coefficients. The direct prime kernel pairs
inputs whose twiddle factors are complex conjugates.
The library requires Rust 1.85 or newer, supports `no_std`, and uses stack
storage without heap allocation. The 3- and 5-point `f32` and 8-point `f32` and `f64`
FFTs use NEON on little-endian AArch64 and SSE2 on x86-64 when enabled for the
target. Both backends share a portable scalar fallback. The crate's unsafe
code is confined to the SIMD module.

## Usage

The top-level functions `fft` and `ifft` accept an array, slice, or vector
of exactly `N` samples and return an array of length `N`.

```
use monarch_butterfly::{fft, ifft};
use num_complex::Complex;

let input: Vec<_> = (0..8).map(|i| Complex::new(i as f32, 0.0)).collect();
let spectrum = fft::<8, _, _>(input);
let scaled_input = ifft::<8, _, _>(spectrum);
```

The inverse is unnormalized. Applying `ifft` to the result of `fft` scales
the original samples by `N`.

Enable the optional `std` Cargo feature to use fused multiply-add in larger
direct prime butterflies when the compilation target supports hardware FMA.

The original FFT kernels were inspired by Elliott Mahler’s [RustFFT](https://github.com/ejmahler/RustFFT).
