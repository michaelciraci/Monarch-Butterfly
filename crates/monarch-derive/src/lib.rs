//! This crate is not meant to be used on its own. This crate containsthe under the hood procedural macros to autogenerate
//! FFTs for [monarch-butterfly](https://crates.io/crates/monarch-butterfly)

use num_complex::Complex;
use num_integer::Integer;
use num_traits::{Float, FloatConst};
use proc_macro::{Span, TokenStream};
use quote::quote;
use syn::Ident;

const SIZES: [usize; 134] = [
    2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30,
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
    55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78,
    79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
    102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120,
    121, 122, 123, 124, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140
];

const HAND_GEN: [usize; 5] = [3, 9, 18, 27, 125];

#[derive(PartialEq, Debug)]
enum FFTType {
    PowerOfTwo,
    Prime,
    Coprime,
    Mixed,
}

impl FFTType {
    fn compute_type(n: usize) -> FFTType {
        if n.is_power_of_two() {
            FFTType::PowerOfTwo
        } else {
            // Check if it has an integer square root
            let sqr = (n as f64).sqrt().round() as usize;
            if sqr * sqr == n && sqr % 2 == 1 {
                FFTType::Mixed
            } else {
                let (v1, _v2) = compute_coprimes(n);
                if v1 == 1 {
                    FFTType::Prime
                } else {
                    FFTType::Coprime
                }
            }
        }
    }
}

fn compute_coprimes(n: usize) -> (usize, usize) {
    let sqr = (n as f32).sqrt().ceil() as usize;
    for v1 in (1..sqr).rev() {
        let v2 = n / v1;
        if v2 * v1 == n {
            if v1.gcd(&v2) == 1 {
                return (usize::min(v1, v2), usize::max(v1, v2));
            }
        }
    }
    todo!()
}

fn compute_twiddle_forward<T: Float + FloatConst>(index: usize, fft_len: usize) -> Complex<T> {
    let constant = T::from(-2.0).unwrap() * T::PI() / T::from(fft_len).unwrap();
    // index * -2PI / fft_len
    let angle = constant * T::from(index).unwrap();

    Complex::new(angle.cos(), angle.sin())
}

#[proc_macro]
pub fn generate_switch(_input: TokenStream) -> TokenStream {
    let mut all_sizes: Vec<_> = SIZES
        .clone()
        .into_iter()
        .chain(HAND_GEN.clone().into_iter())
        .collect();
    all_sizes.sort();

    let ss_forward = all_sizes.clone().into_iter().map(|s| {
        let func = Ident::new(&format!("fft{}", s), Span::call_site().into());

        quote! {
            #s => {
                let x = #func(x_in);
                std::array::from_fn(|i| x[i])
             },
        }
    });
    let ss_inverse = all_sizes.into_iter().map(|s| {
        let func = Ident::new(&format!("ifft{}", s), Span::call_site().into());

        quote! {
            #s => {
                let x = #func(x_in);
                std::array::from_fn(|i| x[i])
             },
        }
    });

    let expanded = quote! {
        #[inline]
        pub fn fft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; N] {
            let x_in = input.as_ref();
            assert_eq!(x_in.len(), N);

            match N {
                1 => { std::array::from_fn(|i| x_in[i]) },
                #(#ss_forward)*
                _ => unimplemented!(),
            }
        }

        #[inline]
        pub fn ifft<const N: usize, T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; N] {
            let x_in = input.as_ref();
            assert_eq!(x_in.len(), N);

            match N {
                1 => { std::array::from_fn(|i| x_in[i]) },
                #(#ss_inverse)*
                _ => unimplemented!(),
            }
        }
    };
    proc_macro::TokenStream::from(expanded)
}

#[proc_macro]
pub fn generate_powers_of_two(_input: TokenStream) -> TokenStream {
    let sizes = SIZES
        .clone()
        .into_iter()
        .filter(|n| FFTType::compute_type(*n) == FFTType::PowerOfTwo);
    let ss = sizes.map(|s| {
        let func = Ident::new(&format!("fft{}", s), Span::call_site().into());
        let half = s / 2;
        let half_butterfly = Ident::new(&format!("fft{}", half), Span::call_site().into());
        let half_butterfly_even_idx = (0..s).step_by(2).map(|f|{
            quote! {
                x[#f],
            }
        });
        let half_butterfly_odd_idx = (1..s).step_by(2).map(|f|{
            quote! {
                x[#f],
            }
        });

        let t_s = (0..half).map(|tt|
            quote! {
                Complex::exp(Complex::<T>::i() * T::from(-2.0).unwrap() * T::PI() * T::from(#tt).unwrap() / T::from(n).unwrap()) * odd[#tt]
            }
        );

        let sum_halves = (0..half).map(|t_e| quote! {
            even[#t_e] + t[#t_e],
        }
        );
        let sub_halves = (0..half).map(|t_o| quote! {
            even[#t_o] - t[#t_o],
        });

        quote! {
            #[inline]
            pub fn #func<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; #s] {
                let n = #s;
                let x = input.as_ref();
                assert_eq!(n, x.len());

                let even: [Complex<T>; #half] = #half_butterfly([
                    #(#half_butterfly_even_idx)*
                ]);
                let odd: [Complex<T>; #half] = #half_butterfly([
                    #(#half_butterfly_odd_idx)*
                ]);

                let t: [Complex<T>; #half] = [
                    #(#t_s),*
                ];

                [
                    #(#sum_halves)*
                    #(#sub_halves)*
                ]
            }
        }
    });

    let expanded = quote! {
        #[inline]
        pub fn fft1<T: Float, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 1] {
            let n = 1;
            let x = input.as_ref();
            assert_eq!(n, x.len());

            [x[0]]
        }

        #(#ss)*
    };
    proc_macro::TokenStream::from(expanded)
}

#[proc_macro]
pub fn generate_coprimes(_input: TokenStream) -> TokenStream {
    let sizes = SIZES
        .clone()
        .into_iter()
        .filter(|n| FFTType::compute_type(*n) == FFTType::Coprime);
    let ss = sizes.map(|s| {
        let (c1, c2) = compute_coprimes(s);
        let func = Ident::new(&format!("fft{}", s), Span::call_site().into());
        let func1 = Ident::new(&format!("fft{}", c1), Span::call_site().into());
        let func2 = Ident::new(&format!("fft{}", c2), Span::call_site().into());

        let rows = (0..c2).map(|i|  {
            let mut start = c1 * i;
            let idx = (0..c1).map(|_| {
                let index = start;
                start = (start + c2) % s;
                quote! {
                    x[#index],
                }}
            );
            let row_call = Ident::new(&format!("row{}", i), Span::call_site().into());

            quote! {
                let #row_call = #func1([ #(#idx)* ]);
        }});

        let cols = (0..c1).map(|i| {
            let idx = (0..c2).map(|ii| {
                let row_call = Ident::new(&format!("row{}", ii), Span::call_site().into());
                quote! {
                    #row_call[#i]
                }
            });

            let col_call = Ident::new(&format!("col{}", i), Span::call_site().into());

            quote! {
                let #col_call = #func2([ #(#idx),*]);
            }
        });

        let combine = (0..s).map(|i| {
            let col = i % c1;
            let idx = i % c2;
            let f = Ident::new(&format!("col{}", col), Span::call_site().into());
            quote! {
                #f[#idx],
            }
        });

        quote! {
            #[inline]
            pub fn #func<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; #s] {
                let n = #s;
                let x = input.as_ref();
                assert_eq!(n, x.len());

                #(#rows)*
                #(#cols)*


                [#(#combine)*]

            }
        }
    });

    let expanded = quote! {
        #(#ss)*
    };
    proc_macro::TokenStream::from(expanded)
}

#[proc_macro]
pub fn generate_mixed_radix(_input: TokenStream) -> TokenStream {
    let sizes = SIZES
        .clone()
        .into_iter()
        .filter(|n| FFTType::compute_type(*n) == FFTType::Mixed);
    let ss = sizes.map(|s| {
        let c1 = (s as f64).sqrt().round() as usize;
        let c2 = c1;
        let func = Ident::new(&format!("fft{}", s), Span::call_site().into());
        let func1 = Ident::new(&format!("fft{}", c1), Span::call_site().into());
        let func2 = Ident::new(&format!("fft{}", c2), Span::call_site().into());

        let rows = (0..c1).map(|i|  {
            let idx = (i..s).step_by(c1).map(|xx| {
                let index = xx % s;
                quote! {
                    x[#index],
                }}
            );
            let row_call = Ident::new(&format!("row{}", i), Span::call_site().into());

            quote! {
                let #row_call = #func2([ #(#idx)* ]);
        }});

        let mut twiddles = vec![Complex::<f64>::new(0.0, 0.0); s];
        for (x, twiddle_chunk) in twiddles.chunks_exact_mut(c2).enumerate() {
            for (y, twiddle_element) in twiddle_chunk.iter_mut().enumerate() {
                *twiddle_element = compute_twiddle_forward(x * y, s);
            }
        }

        let cols = (0..c2).map(|i| {
            let mut start_idx = i;
            let idx = (0..c1).map(|ii| {
                let row_call = Ident::new(&format!("row{}", ii), Span::call_site().into());
                let re = twiddles[start_idx].re;
                let im = twiddles[start_idx].im;
                start_idx += c2;
                quote! {
                    #row_call[#i] * Complex::new(T::from(#re).unwrap(), T::from(#im).unwrap())
                }
            });

            let col_call = Ident::new(&format!("col{}", i), Span::call_site().into());

            quote! {
                let #col_call = #func1([ #(#idx),*]);
            }
        });

        let combine = (0..s).map(|i| {
            let col = i % c2;
            let idx = i / c2;
            let f = Ident::new(&format!("col{}", col), Span::call_site().into());
            quote! {
                #f[#idx],
            }
        });

        quote! {
            #[inline]
            pub fn #func<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; #s] {
                let n = #s;
                let x = input.as_ref();
                assert_eq!(n, x.len());

                #(#rows)*
                #(#cols)*


                [#(#combine)*]

            }
        }
    });

    let expanded = quote! {
        #(#ss)*
    };
    proc_macro::TokenStream::from(expanded)
}

#[proc_macro]
pub fn generate_primes(_input: TokenStream) -> TokenStream {
    let sizes = SIZES
        .clone()
        .into_iter()
        .filter(|n| FFTType::compute_type(*n) == FFTType::Prime);
    let ss = sizes.map(|s| {
        let func = Ident::new(&format!("fft{}", s), Span::call_site().into());
        let halflen = (s + 1) / 2;
        let twiddles = (1..halflen).map(|n| {
            let var = Ident::new(&format!("twiddle{}", n), Span::call_site().into());
            let val: Complex<f64> = compute_twiddle_forward(n, s);
            let re = val.re;
            let im = val.im;
            quote! {
                let #var = Complex::new(T::from(#re).unwrap(), T::from(#im).unwrap());
            }
        });
        let first_codegen = (1..halflen).map(|n| {
            let var1 = Ident::new(&format!("x{}{}p", n, s - n), Span::call_site().into());
            let var2 = Ident::new(&format!("x{}{}n", n, s - n), Span::call_site().into());
            quote! {
                let #var1 = x[#n] + x[#s - #n];
                let #var2 = x[#n] - x[#s - #n];
            }
        });
        let second_codegen = (1..halflen).map(|n| {
            let var = Ident::new(&format!("x{}{}p", n, s - n), Span::call_site().into());
            quote! {
                + #var
            }
        });
        let third_codegen = (1..halflen).map(|n| {
            let var1 = Ident::new(&format!("b{}{}re_a", n, s - n), Span::call_site().into());
            let sub1 = (1..halflen).map(|m| {

                let mut mn = (m * n) % s;
                if mn > s / 2 {
                    mn = s - mn;
                }
                let var2 = Ident::new(&format!("twiddle{}", mn), Span::call_site().into());
                let var3 = Ident::new(&format!("x{}{}p", m, s - m), Span::call_site().into());
                quote! {
                    + #var2.re * #var3.re
                }
            });

            quote! {
                let #var1 = x[0].re #(#sub1)* ;
            }
        });

        let fourth_codegen = (1..halflen).map(|n| {
            let var1 = Ident::new(&format!("b{}{}re_b", n, s - n), Span::call_site().into());
            let sub1 = (1..halflen).map(|m| {

                let mut mn = (m * n) % s;
                if mn > s / 2 {
                    mn = s - mn;
                    let var2 = Ident::new(&format!("twiddle{}", mn), Span::call_site().into());
                    let var3 = Ident::new(&format!("x{}{}n", m, s - m), Span::call_site().into());
                    quote! {
                        - #var2.im * #var3.im
                    }
                } else {
                    let var2 = Ident::new(&format!("twiddle{}", mn), Span::call_site().into());
                    let var3 = Ident::new(&format!("x{}{}n", m, s - m), Span::call_site().into());
                    quote! {
                        + #var2.im * #var3.im
                    }
                }
            });

            quote! {
                let #var1 = T::zero() #(#sub1)* ;
            }
        });
        let fifth_codegen = (1..halflen).map(|n| {
            let var1 = Ident::new(&format!("b{}{}im_a", n, s - n), Span::call_site().into());
            let sub1 = (1..halflen).map(|m| {

                let mut mn = (m * n) % s;
                if mn > s / 2 {
                    mn = s - mn;
                }
                let var2 = Ident::new(&format!("twiddle{}", mn), Span::call_site().into());
                let var3 = Ident::new(&format!("x{}{}p", m, s - m), Span::call_site().into());
                quote! {
                    + #var2.re * #var3.im
                }
            });

            quote! {
                let #var1 = x[0].im #(#sub1)* ;
            }
        });

        let sixth_codegen = (1..halflen).map(|n| {
            let var1 = Ident::new(&format!("b{}{}im_b", n, s - n), Span::call_site().into());
            let sub1 = (1..halflen).map(|m| {

                let mut mn = (m * n) % s;
                if mn > s / 2 {
                    mn = s - mn;
                    let var2 = Ident::new(&format!("twiddle{}", mn), Span::call_site().into());
                    let var3 = Ident::new(&format!("x{}{}n", m, s - m), Span::call_site().into());
                    quote! {
                        - #var2.im * #var3.re
                    }
                } else {
                    let var2 = Ident::new(&format!("twiddle{}", mn), Span::call_site().into());
                    let var3 = Ident::new(&format!("x{}{}n", m, s - m), Span::call_site().into());
                    quote! {
                        + #var2.im * #var3.re
                    }
                }
            });

            quote! {
                let #var1 = T::zero() #(#sub1)* ;
            }
        });

        let seventh_codegen = (1..s).map(|n| {
            let mut nfold = n;
            if n > s / 2 {
                nfold = s - n;
                let var1 = Ident::new(&format!("out{}re", n), Span::call_site().into());
                let var2 = Ident::new(&format!("out{}im", n), Span::call_site().into());
                let var3 = Ident::new(&format!("b{}{}re_a", nfold, s-nfold), Span::call_site().into());
                let var4 = Ident::new(&format!("b{}{}re_b", nfold, s-nfold), Span::call_site().into());
                let var5 = Ident::new(&format!("b{}{}im_a", nfold, s-nfold), Span::call_site().into());
                let var6 = Ident::new(&format!("b{}{}im_b", nfold, s-nfold), Span::call_site().into());
                quote! {
                    let #var1 = #var3 + #var4;
                    let #var2 = #var5 - #var6;
                }
            } else {
                let var1 = Ident::new(&format!("out{}re", n), Span::call_site().into());
                let var2 = Ident::new(&format!("out{}im", n), Span::call_site().into());
                let var3 = Ident::new(&format!("b{}{}re_a", nfold, s-nfold), Span::call_site().into());
                let var4 = Ident::new(&format!("b{}{}re_b", nfold, s-nfold), Span::call_site().into());
                let var5 = Ident::new(&format!("b{}{}im_a", nfold, s-nfold), Span::call_site().into());
                let var6 = Ident::new(&format!("b{}{}im_b", nfold, s-nfold), Span::call_site().into());
                quote!{
                    let #var1 = #var3 - #var4;
                    let #var2 = #var5 + #var6;
                }
            }
        });
        let eigth_codegen = (1..s).map(|n| {
            let var_re: Ident = Ident::new(&format!("out{}re", n), Span::call_site().into());
            let var_im: Ident = Ident::new(&format!("out{}im", n), Span::call_site().into());
            quote! {
                Complex::new(#var_re, #var_im),
            }
        });

        quote! {
            #[inline]
            pub fn #func<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; #s] {
                let n = #s;
                let x = input.as_ref();
                assert_eq!(n, x.len());

                #(#twiddles)*

                #(#first_codegen)*
                let sum = x[0] #(#second_codegen)* ;
                #(#third_codegen)*
                #(#fourth_codegen)*
                #(#fifth_codegen)*
                #(#sixth_codegen)*
                #(#seventh_codegen)*

                [
                    sum,
                    #(#eigth_codegen)*
                ]
            }

    }});

    let expanded = quote! {
        #(#ss)*
    };
    proc_macro::TokenStream::from(expanded)
}

#[proc_macro]
pub fn generate_iffts(_input: TokenStream) -> TokenStream {
    let mut all_sizes: Vec<_> = SIZES
        .clone()
        .into_iter()
        .chain(HAND_GEN.clone().into_iter())
        .collect();
    all_sizes.sort();
    let iffts = all_sizes.into_iter().map(|n| {
        let func = Ident::new(&format!("ifft{}", n), Span::call_site().into());
        let input_args = (0..n).map(|i| {

            quote! {
                x[#i].conj(),
            }
        });
        let output_args = (0..n).map(|i| {

            quote! {
                out[#i].conj(),
            }
        });

        quote! {
            #[inline]
            pub fn #func<T: Float + FloatConst, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; #n] {
                let x = input.as_ref();
                assert_eq!(x.len(), #n);

                let out: [Complex<T>; #n] = fft::<#n, _, _>([
                    #(#input_args)*
                ]);
                [
                    #(#output_args)*
                ]

            }
        }
    });

    let expanded = quote! {
        #[inline]
        pub fn ifft1<T: Float, A: AsRef<[Complex<T>]>>(input: A) -> [Complex<T>; 1] {
            let n = 1;
            let x = input.as_ref();
            assert_eq!(n, x.len());

            [x[0]]
        }

        #(#iffts)*
    };
    proc_macro::TokenStream::from(expanded)
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_coprimes() {
        let coprimes = vec![
            (2, 3),
            (2, 5),
            (3, 4),
            (2, 7),
            (3, 5),
            (4, 5),
            (3, 7),
            (2, 11),
            (3, 8),
            (2, 13),
            (4, 7),
            (5, 6),
            (3, 11),
            (2, 17),
            (5, 7),
            (2, 19),
            (3, 13),
            (5, 8),
            (6, 7),
            (4, 11),
            (5, 9),
            (2, 23),
            (2, 25),
            (3, 17),
            (4, 13),
            (2, 27),
            (5, 11),
            (7, 8),
            (3, 19),
            (2, 29),
            (5, 12),
            (2, 31),
            (7, 9),
            (5, 13),
            (6, 11),
            (4, 17),
            (3, 23),
            (7, 10),
            (8, 9),
            (2, 37),
            (3, 25),
            (4, 19),
            (7, 11),
            (6, 13),
            (5, 16),
            (2, 41),
            (7, 12),
            (5, 17),
            (2, 43),
            (3, 29),
        ];
        for (v1, v2) in coprimes {
            let n = v1 * v2;
            let (computed_v1, computed_v2) = compute_coprimes(n);
            dbg!(n, v1, v2, computed_v1, computed_v2);
            assert_eq!(v1, computed_v1);
            assert_eq!(v2, computed_v2);
        }
    }

    #[test]
    fn test_fft_type() {
        assert_eq!(FFTType::compute_type(25), FFTType::Mixed);
        assert_eq!(FFTType::compute_type(36), FFTType::Coprime);
    }
}
