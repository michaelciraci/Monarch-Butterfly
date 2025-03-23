//! This crate is not meant to be used on its own. This crate containsthe under the hood procedural macros to autogenerate
//! FFTs for [monarch-butterfly](https://crates.io/crates/monarch-butterfly)

use proc_macro::{Span, TokenStream};
use quote::quote;
use syn::Ident;

const SIZES: [usize; 11] = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];

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
