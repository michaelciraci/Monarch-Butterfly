[![Build](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml/badge.svg)](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml)
[![](https://img.shields.io/crates/v/monarch-butterfly)](https://crates.io/crates/monarch-butterfly)
[![](https://docs.rs/monarch-butterfly/badge.svg)](https://docs.rs/monarch-butterfly/)

# Monarch Butterfly

This prototype computes forward and inverse FFTs for 1 through 10 complex
samples. The size is a const generic, and twiddle factors are computed at
compile time on stable Rust. A size dispatcher selects fixed-size kernels or
composes smaller transforms. The composition takes its factors as separate
const parameters, so intermediate arrays have exact lengths without unstable
generic const expressions. The 5-point kernel shares work between paired bins;
the 7-point leaf pairs conjugate terms in a direct DFT.

The library supports `no_std` and uses stack storage without heap allocation.
It is an example implementation, not a performance-optimized release.

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
