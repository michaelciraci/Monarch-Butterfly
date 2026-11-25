use std::{env, fmt::Write, fs, path::PathBuf};

// Keep code size bounded. The library handles larger lengths iteratively.
const MAX_KERNEL_LEN: usize = 200;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let mut kernels = String::new();
    let mut dispatch = String::new();
    for n in 2..=MAX_KERNEL_LEN {
        let mut rest = n >> n.trailing_zeros();
        while rest % 3 == 0 {
            rest /= 3;
        }
        if rest == 1 {
            generate_kernel(&mut kernels, n);
            writeln!(dispatch, "{n} => out.copy_from_slice(&fft_{n}(input)),").unwrap();
        }
    }
    writeln!(
        kernels,
        "
        const MAX_KERNEL_LEN: usize = {MAX_KERNEL_LEN};
        #[inline(always)]
        fn small_fft<const N: usize, T: Float + FloatConst + ConstZero>(
            input: &[Complex<T>],
        ) -> [Complex<T>; N] {{
            let mut out = [Complex::<T>::ZERO; N];
            match N {{ {dispatch} _ => panic!(\"unsupported FFT length {{}}\", N) }}
            out
        }}"
    )
    .unwrap();
    let dir = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    fs::write(dir.join("kernels.rs"), kernels).unwrap();
}

fn generate_kernel(code: &mut String, n: usize) {
    writeln!(
        code,
        "#[inline(always)]
        #[allow(clippy::identity_op)] // Generated offsets.
        fn fft_{n}<T: Float + FloatConst + ConstZero>(x: &[Complex<T>]) -> [Complex<T>; {n}] {{"
    )
    .unwrap();
    if n <= 4 {
        let args = indices(n, |i| format!("x[{i}]"));
        writeln!(code, "radix{n}_butterfly({args})\n}}").unwrap();
        return;
    }

    // Factor into two small transforms of similar size. Coprime factors
    // allow a Good-Thomas permutation to eliminate the inter-stage twiddles.
    let cols = (2..)
        .take_while(|c| c * c <= n)
        .filter(|c| n % c == 0)
        .last()
        .unwrap();
    let rows = n / cols;
    // The measured faster layouts use coprime factors at 24 and four
    // independent eight-point transforms at 32.
    let (rows, cols) = match n {
        24 => (3, 8),
        32 => (4, 8),
        _ => (rows, cols),
    };
    let coprime = (2..=cols).all(|d| rows % d != 0 || cols % d != 0);
    let inputs = indices(cols, |i| {
        let stride = i * rows;
        if coprime {
            format!("x[(r * {cols} + {stride}) % {n}]")
        } else {
            format!("x[r + {stride}]")
        }
    });
    writeln!(
        code,
        "let mut sub = [[Complex::<T>::ZERO; {rows}]; {cols}];
        for r in 0..{rows} {{
            let column = fft_{cols}(&[{inputs}]);
            for k in 0..{cols} {{ sub[k][r] = column[k]; }}
        }}"
    )
    .unwrap();

    for k in 0..cols {
        let inputs = indices(rows, |r| {
            let value = format!("sub[{k}][{r}]");
            if coprime {
                value
            } else {
                twisted(value, k * r, n)
            }
        });
        writeln!(code, "let row{k} = fft_{rows}(&[{inputs}]);").unwrap();
    }

    let outputs = indices(n, |i| {
        let row = i % cols;
        let bin = if coprime { i % rows } else { i / cols };
        format!("row{row}[{bin}]")
    });
    writeln!(code, "[{outputs}]\n}}").unwrap();
}

fn indices(count: usize, f: impl FnMut(usize) -> String) -> String {
    (0..count).map(f).collect::<Vec<_>>().join(", ")
}

fn twisted(value: String, power: usize, n: usize) -> String {
    match (4 * power / n, 4 * power % n) {
        (0, 0) => value,
        (1, 0) => format!("Complex::new({value}.im, -{value}.re)"),
        (2, 0) => format!("-{value}"),
        (3, 0) => format!("Complex::new(-{value}.im, {value}.re)"),
        _ => {
            let w = format!("cast_twiddle::<T>(const {{ forward_twiddle({power}, {n}) }})");
            if n <= 8 {
                format!("{value} * {w}")
            } else {
                format!("{{ let w = {w}; rotate({value}, w.re, w.im) }}")
            }
        }
    }
}
