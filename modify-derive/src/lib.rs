use proc_macro::TokenStream;
use quote::{format_ident, quote};
use std::iter::Iterator;
use syn::{parse_macro_input, DeriveInput, Fields, Ident, Type};

#[proc_macro_derive(Modify)]
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
                    let syn::Field { ident: id, ty, .. } = nf;
                    fids.push(id.unwrap());
                    ftys.push(ty);
                }
            }
            _ => panic!("overrides only apply to named fields"),
        },
        _ => panic!(
            "Override macro expects struct, but found enum/union"
        ),
    };

    let output = quote!(
        impl #id {
            #(pub fn modify_#fids(mut self, new_value: #ftys) -> #id {
                self.#fids = new_value;
                self
            })*
        }
    );

    TokenStream::from(output)
}
