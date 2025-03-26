## Monarch Butterfly

Experimental FFT library with an emphasis on runtime performance.
The goal of this library is to hand-unroll all loops in a procedural macro
for optimal SIMD throughput.

This currently only works on powers of two.

This library implements FFTs for both `f32` and `f64` for the following sizes:
```
1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048
```

This library will use all SIMD features your CPU has available including AVX512,
assuming you compile with those features (`RUSTFLAGS="-C target-cpu=native" cargo build`).

The larger the FFT sizes, the larger speed boost this library will give you.

As an example of AVX512 instructions, here is an example on just an FFT
of size 128: https://godbolt.org/z/rz48azEsd (`Ctrl+F` for "zmm" instructions)

If a larger FFT size is needed, just clone the repo and add the needed
sizes to the top of `crates\monarch-derive\src\lib.rs` and larger FFTs
will be generated. However, this comes at the cost of a longer compile time.

```
use monarch_butterfly::*;
use num_complex::Complex;

let input: Vec<_> = (0..8).map(|i| Complex::new(i as f32, 0.0)).collect();
let output_slice = fft8(&input);
let output_vec = fft8(input);
```

The FFTs before unrolling are heavily inspired from [RustFFT](https://github.com/ejmahler/RustFFT).
Credit is given to Elliott Mahler as the RustFFT original author.
