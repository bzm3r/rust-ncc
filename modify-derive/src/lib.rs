// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use proc_macro2::Ident;
use quote::quote;
use syn::DeriveInput;
use syn::{parse_macro_input, Type};

#[proc_macro_derive(Modify)]
pub fn derive(
    input: proc_macro::TokenStream,
) -> proc_macro::TokenStream {
    let DeriveInput {
        ident: id, data, ..
    } = parse_macro_input!(input as DeriveInput);

    let mut fids: Vec<Ident> = vec![];
    let mut ftys: Vec<Type> = vec![];
    match data {
        syn::Data::Struct(syn::DataStruct {
            fields: syn::Fields::Named(syn::FieldsNamed { named, .. }),
            ..
        }) => {
            for nf in named {
                let syn::Field { ident: id, ty, .. } = nf;
                fids.push(id.unwrap());
                ftys.push(ty);
            }
        },
        _ => panic!("can only derive modify methods for named fields of structs"),
    };

    let method_names = fids.iter().map(|fid| {
        let m_name = format!("modify_{}", fid);
        syn::Ident::new(&m_name, fid.span())
    });

    let methods = method_names.zip(fids.iter().zip(ftys)).map(
        |(mn, (fid, fty))| {
            quote!(
                pub fn #mn(mut self, val: #fty) -> Self {
                    self.#fid = val;
                    self
                }
            )
        },
    );

    let expanded = quote! {
        impl #id {
            #(#methods)*
        }
    };

    expanded.into()
}
