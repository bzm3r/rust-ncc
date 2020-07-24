mod avro_schema;
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

    let fid_strings: Vec<String> = fids.iter().map(|fid| fid.to_string()).collect();
    let fschemas: Vec<proc_macro2::TokenStream> = ftys
        .iter().zip(fid_strings.iter())
        .map(|(fty, fstr)| {
            avro_schema::from_syn(fstr.as_str(), fty)
        })
        .collect();
    let fpositions = fids.iter().enumerate().map(|(i, _)| i).collect::<Vec<usize>>();

    proc_macro::TokenStream::from(quote!(
        impl #id {
            pub fn schematize(prev_namespace: Option<String>) -> avro_rs::schema::Schema {
                let new_namespace = if let Some(pns) = prev_namespace.as_ref() {
                    Some(vec![String::from(pns), String::from(#id_string)].join("."))
                } else {
                    Some(String::from(#id_string))
                };
                avro_rs::schema::Schema::Record {
                    name: avro_rs::schema::Name {
                        name: String::from(#id_string),
                        namespace: prev_namespace.clone(),
                        aliases: None,
                    },
                    doc: None,
                    fields: vec![#(
                        avro_rs::schema::RecordField {
                            name: std::string::String::from(#fid_strings),
                            doc: None,
                            default: None,
                            schema: #fschemas,
                            order: avro_rs::schema::RecordFieldOrder::Ignore,
                            position: #fpositions,
                        }
                    ),*],
                    lookup: {
                        let mut r: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
                        #(r.insert(std::string::String::from(#fid_strings), #fpositions);)*
                        r
                    },
                }
            }
        }
    ))
}
