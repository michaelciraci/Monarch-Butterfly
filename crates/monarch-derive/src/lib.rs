//! This crate is not meant to be used on its own. This crate containsthe under the hood procedural macros to autogenerate
//! FFTs for [monarch-butterfly](https://crates.io/crates/monarch-butterfly)

use num_complex::Complex;
use num_traits::{Float, FloatConst};
use proc_macro::{Span, TokenStream};
use quote::quote;
use syn::Ident;

const SIZES: [usize; 11] = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
const COPRIMES: [(usize, usize); 11] = [
    (2, 5),
    (4, 3),
    (2, 7),
    (3, 5),
    (4, 5),
    (3, 7),
    (2, 11),
    (3, 8),
    (2, 13),
    (4, 7),
    (5, 6),
];
const MIXED_RADIX: [(usize, usize); 1] = [(5, 5)];
const PRIMES: [usize; 9] = [5, 7, 11, 13, 17, 19, 23, 29, 31];

fn _compute_twiddle<T: Float + FloatConst>(index: usize, fft_len: usize) -> Complex<T> {
    let constant = T::from(-2.0).unwrap() * T::PI() / T::from(fft_len).unwrap();
    // index * -2PI / fft_len
    let angle = constant * T::from(index).unwrap();

    Complex::new(angle.cos(), angle.sin())
}

#[proc_macro]
pub fn generate_powers_of_two(_input: TokenStream) -> TokenStream {
    let ss = SIZES.clone().into_iter().map(|s| {
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
        pub fn fft1<T: Float>(x: [Complex<T>; 1]) -> [Complex<T>; 1] {
            x
        }

        #(#ss)*
    };
    proc_macro::TokenStream::from(expanded)
}

#[proc_macro]
pub fn generate_coprimes(_input: TokenStream) -> TokenStream {
    let ss = COPRIMES.clone().into_iter().map(|(c1, c2)| {
        let s = c1 * c2;
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
    let ss = MIXED_RADIX.clone().into_iter().map(|(c1, c2)| {
        let s = c1 * c2;
        let func = Ident::new(&format!("fft{}", s), Span::call_site().into());
        let func1 = Ident::new(&format!("fft{}", c1), Span::call_site().into());
        let func2 = Ident::new(&format!("fft{}", c2), Span::call_site().into());

        let rows = (0..c2).map(|i|  {
            // let mut start = c1 * i;
            let idx = (i..s).step_by(c1).map(|xx| {
                let index = xx % s;
                // start = (start + c2) % s;
                quote! { 
                    x[#index],  
                }}
            );
            let row_call = Ident::new(&format!("row{}", i), Span::call_site().into());
            
            quote! {
                let #row_call = #func1([ #(#idx)* ]);
        }});

        let mut twiddles = vec![Complex::<f64>::new(0.0, 0.0); 25];
        for (x, twiddle_chunk) in twiddles.chunks_exact_mut(5).enumerate() {
            for (y, twiddle_element) in twiddle_chunk.iter_mut().enumerate() {
                *twiddle_element = _compute_twiddle(x * y, 25);
            }
        }

        let cols = (0..c1).map(|i| {
            let mut start_idx = i;
            let idx = (0..c2).map(|ii| {
                let row_call = Ident::new(&format!("row{}", ii), Span::call_site().into());
                let re = twiddles[start_idx].re;
                let im = twiddles[start_idx].im;
                start_idx = (start_idx + c1) % s;
                quote! {
                    #row_call[#i] * Complex::new(T::from(#re).unwrap(), T::from(#im).unwrap())
                }
            });

            let col_call = Ident::new(&format!("col{}", i), Span::call_site().into());

            quote! {
                let #col_call = #func2([ #(#idx),*]);
            }
        });

        let combine = (0..s).map(|i| {
            let col = i % c1;
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
    let ss = PRIMES.clone().into_iter().map(|s| {
        let func = Ident::new(&format!("fft{}", s), Span::call_site().into());
        let halflen = (s + 1) / 2;
        let twiddles = (1..halflen).map(|n| {
            let var = Ident::new(&format!("twiddle{}", n), Span::call_site().into());
            let val: Complex<f64> = _compute_twiddle(n, s);
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
