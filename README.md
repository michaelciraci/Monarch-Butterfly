[![Build](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml/badge.svg)](https://github.com/michaelciraci/Monarch-Butterfly/actions/workflows/rust.yml)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success.svg)](https://github.com/rust-secure-code/safety-dance/)
[![](https://img.shields.io/crates/v/monarch-butterfly)](https://crates.io/crates/monarch-butterfly)
[![](https://docs.rs/monarch-butterfly/badge.svg)](https://docs.rs/monarch-butterfly/)

# Monarch Butterfly

Experimental FFT library where all FFTs are proc-macro generated const-evaluation functions. 
The one requirement is you must know the size of the FFT at compile-time. Knowing the FFT size at
compile time gives immense gains, as the compiler is able unroll the call stack and optimize for
SIMD throughput through function calls.

This library implements FFTs for both `f32` and `f64` sizes `1-200`. The FFTs are auto-generated so this limit could be increased above 200 at the expense of compile time.

## Features

- All functions are auto-generated with proc-macros with unrolled loops
- Zero `unsafe` code
- no_std
- Completely portable
- Const-evaluation functions

## Limitations

- FFT size must be known at compile time
- By default, only FFTs up to size 200 are generated

## Comparison

RustFFT and FFTW FFT sizes can be decided at runtime, so comparing against a library whose FFT sizes need to be known at compile time is comparing apples and oranges. With that in mind, knowing the FFT size at compile time does give immense gains.

![log](assets/log_comparison.png)

## Usage

The top level functions are `fft` and `ifft`.

```
use monarch_butterfly::*;
use num_complex::Complex;

let input: Vec<_> = (0..8).map(|i| Complex::new(i as f32, 0.0)).collect();
let output = fft::<8, _, _>(input);
```

This library will use all SIMD features your CPU has, assuming `rustc` can compile to those SIMD features.

The larger the FFT sizes, the larger speed boost this library will give you.

As an example of AVX512 instructions, here is an example on just an FFT
of size 128: <https://godbolt.org/z/Y58eh1x5a> (`Ctrl+F` for "zmm" instructions)

The FFTs before unrolling are heavily inspired from [RustFFT](https://github.com/ejmahler/RustFFT).
Credit is given to Elliott Mahler as the RustFFT original author.
