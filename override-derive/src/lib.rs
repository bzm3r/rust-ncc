use proc_macro::TokenStream;
use quote::{format_ident, quote};
use std::iter::Iterator;
use syn::{
    parse_macro_input, DeriveInput, Fields, Ident, Type,
};

#[proc_macro_derive(Overrides)]
pub fn derive(input: TokenStream) -> TokenStream {
    let DeriveInput {
        ident: id, data, ..
    } = parse_macro_input!(input as DeriveInput);

    let mut fids: Vec<Ident> = vec![];
    let mut ftys: Vec<Type> = vec![];
    match data {
        syn::Data::Struct(ds) => match ds.fields {
            Fields::Named(nfs) => {
                for nf in nfs.named {
                    let syn::Field {
                        ident: id,
                        ty,
                        ..
                    } = nf;
                    fids.push(id.unwrap());
                    ftys.push(ty);
                }
            }
            _ => panic!("overrides only apply to named fields"),
        },
        _ => panic!("Override macro expects struct, but found enum/union"),
    };

    let oid = format_ident!("{}Overrides", &id.to_string());

    let output = quote!(
        #[derive(serde::Deserialize)]
        pub struct #oid {
            #(#fids: std::option::Option<#ftys>,)*
        }

        impl #id {
            pub fn apply_overrides(&self, overrides: &#oid) -> #id {
                let mut r = Self::default();

                #(if let Some(dat) = &overrides.#fids {
                    r.#fids = dat.clone();
                } else {
                    r.#fids = (&self.#fids).clone();
                })*

                r
            }
        }
    );

    TokenStream::from(output)
}
