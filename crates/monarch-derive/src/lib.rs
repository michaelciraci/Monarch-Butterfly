use proc_macro::{Span, TokenStream};
use quote::quote;
use syn::Ident;

const SIZES: [usize; 13] = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];

#[proc_macro]
pub fn generate_powers_of_two(_input: TokenStream) -> TokenStream {
    let ss = SIZES.clone().into_iter().map(|s| {
        let func = Ident::new(&format!("butterfly{}", s), Span::call_site().into());
        let half = s / 2;
        let half_butterfly = Ident::new(&format!("butterfly{}", half), Span::call_site().into());
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
                Complex::exp(-2.0 * Complex::<f32>::i() * std::f32::consts::PI * #tt as f32 / n as f32) * odd[#tt]
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
            pub fn #func(x: [Complex<f32>; #s]) -> [Complex<f32>; #s] {
                let n = #s;

                let even: [Complex<f32>; #half] = #half_butterfly([
                    #(#half_butterfly_even_idx)*
                ]);
                let odd: [Complex<f32>; #half] = #half_butterfly([
                    #(#half_butterfly_odd_idx)*
                ]);

                let t: [Complex<f32>; #half] = [
                    #(#t_s),*
                ];

                [
                    #(#sum_halves)*
                    #(#sub_halves)*
                ]
            }
        }
    });
    // let s = 2;
    // let func: Ident = format!("butterfly{}", s).parse().unwrap();

    let expanded = quote! {
        #[inline]
        pub fn butterfly1(x: [Complex<f32>; 1]) -> [Complex<f32>; 1] {
            x
        }

        #(#ss)*
    };
    proc_macro::TokenStream::from(expanded)
}
