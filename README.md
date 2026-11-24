[![Build](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml/badge.svg)](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success.svg)](https://github.com/rust-secure-code/safety-dance/)
[![](https://img.shields.io/crates/v/monarch-butterfly)](https://crates.io/crates/monarch-butterfly)
[![](https://docs.rs/monarch-butterfly/badge.svg)](https://docs.rs/monarch-butterfly/)

# Monarch Butterfly

Monarch Butterfly computes forward and inverse FFTs of 1 through 200 complex
samples in `f32` or `f64`. The transform size is a const generic. Power-of-two
sizes use iterative radix-2 stages with compile-time butterfly widths. Other
sizes use mixed-radix stages whose lengths, radices, and strides are derived
from `N` at compile time. The FFT twiddle factors are also computed at compile
time. Prime stages pair inputs whose twiddle factors are complex conjugates.
The library requires Rust 1.85 or newer, supports `no_std`, and uses stack
storage without heap allocation or unsafe code.

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

The original FFT kernels were inspired by Elliott Mahler’s [RustFFT](https://github.com/ejmahler/RustFFT).
