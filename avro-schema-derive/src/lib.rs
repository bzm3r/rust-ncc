mod avro_schema;
mod print;

use crate::avro_schema::AvroType;
use quote::quote;
use std::iter::Iterator;
use syn::{parse_macro_input, DeriveInput, Fields, Ident, Type};

#[proc_macro_derive(Schematize)]
pub fn derive(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let DeriveInput {
        ident: id, data, ..
    } = parse_macro_input!(input as DeriveInput);
    let id_string: String = id.to_string();

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
        _ => panic!("Override macro expects struct, but found enum/union"),
    };

    let mut fid_strings: Vec<String> = fids.iter().map(|fid| fid.to_string()).collect();
    let mut fschemas: Vec<proc_macro2::TokenStream> = ftys
        .iter()
        .map(|fty| {
            let pt = avro_schema::RustType::from_syn_type(fty);
            AvroType::from_rust_type(&pt).schema_macro_entry()
        })
        .collect();
    let last_fid_string = fid_strings.remove(fid_strings.len() - 1);
    let last_fschema = fschemas.remove(fschemas.len() - 1);
    let output = quote!(
        impl #id {
            pub fn raw_schema() -> String {
                use std::fmt::Write;
                let mut s = std::string::String::new();

                write!(s, "{{ \"type\": \"record\", \"name\": \"{}\", \"fields\": [ ", #id_string).unwrap();
                #(write!(s, "{{ \"name\": \"{}\", \"type\": {} }}, ", #fid_strings, #fschemas).unwrap();)*
                write!(s, "{{ \"name\": \"{}\", \"type\": {} }} ", #last_fid_string, #last_fschema).unwrap();
                write!(s, "] }}").unwrap();

                return s
            }
        }
    );

    proc_macro::TokenStream::from(output)
}
