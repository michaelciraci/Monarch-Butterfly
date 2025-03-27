//! This crate is not meant to be used on its own. This crate containsthe under the hood procedural macros to autogenerate
//! FFTs for [monarch-butterfly](https://crates.io/crates/monarch-butterfly)

use proc_macro::{Span, TokenStream};
use quote::quote;
use syn::Ident;

const SIZES: [usize; 11] = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
const COPRIMES: [(usize, usize); 9] = [
    (2, 5),
    (4, 3),
    (2, 7),
    (3, 5),
    (4, 5),
    (3, 7),
    (2, 11),
    (3, 8),
    (2, 13),
];

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
