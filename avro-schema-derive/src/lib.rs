mod avro_schema;
use avro_rs::Schema;
use quote::quote;
use std::iter::Iterator;
use syn::{
    parse_macro_input, Attribute, Data, DataEnum, DataStruct,
    DeriveInput, Fields, FieldsNamed, Generics, Ident, Meta,
    NestedMeta, Type, Variant,
};

fn skip_field(attrs: &[Attribute]) -> bool {
    attrs.iter().any(|attr| match attr.parse_meta() {
        Ok(meta) => match meta {
            Meta::List(ml) => {
                ml.path.is_ident("serde")
                    && ml.nested.iter().any(|nm| match nm {
                        NestedMeta::Meta(inner) => {
                            inner.path().is_ident("skip")
                        }
                        _ => false,
                    })
            }
            _ => false,
        },
        _ => false,
    })
}

fn gen_struct_impl(
    id: Ident,
    generics: Generics,
    nfs: FieldsNamed,
) -> proc_macro2::TokenStream {
    let (impl_generics, ty_generics, where_clause) =
        generics.split_for_impl();
    let id_string: String = format!("{}", id.to_string());

    let mut fids: Vec<Ident> = vec![];
    let mut ftys: Vec<Type> = vec![];
    let mut fid_strings: Vec<String> = vec![];
    let mut fschemas: Vec<proc_macro2::TokenStream> = vec![];
    let mut fpositions: Vec<usize> = vec![];

    for (pos, nf) in nfs.named.iter().enumerate() {
        let syn::Field {
            ident, ty, attrs, ..
        } = nf;
        if !skip_field(attrs) {
            let id = ident.as_ref().unwrap().clone();
            let id_str = id.to_string();
            let schema = avro_schema::from_syn(&id_str, ty);
            fids.push(id);
            ftys.push(ty.clone());
            fid_strings.push(id_str);
            fschemas.push(schema);
            fpositions.push(pos);
        }
    }

    quote!(
        impl#impl_generics crate::utils::avro::Schematized for #id #ty_generics #where_clause {
            // Hello world! #impl_generics
            fn schematize(parent_space: Option<String>) -> avro_rs::schema::Schema {
                let new_namespace = if let Some(pns) = parent_space.as_ref() {
                    Some(vec![String::from(pns), String::from(#id_string)].join("."))
                } else {
                    Some(String::from(#id_string))
                };
                avro_rs::schema::Schema::Record {
                    name: avro_rs::schema::Name {
                        name: String::from(#id_string),
                        namespace: parent_space.clone(),
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
    )
}

fn gen_enum_impl(
    id: Ident,
    generics: Generics,
    nfs: Vec<Variant>,
) -> proc_macro2::TokenStream {
    unimplemented!()
    // let (impl_generics, ty_generics, where_clause) =
    //     generics.split_for_impl();
    // let id_string: String = format!("{}", id.to_string());
    //
    // let mut fids: Vec<Ident> = vec![];
    // let mut ftys: Vec<Type> = vec![];
    // let mut fid_strings: Vec<String> = vec![];
    // let mut fschemas: Vec<proc_macro2::TokenStream> = vec![];
    // let mut fpositions: Vec<usize> = vec![];
    //
    // for (pos, nf) in nfs.named.iter().enumerate() {
    //     let syn::Field {
    //         ident: id,
    //         ty,
    //         attrs,
    //         ..
    //     } = nf;
    //     if !skip_field(&attrs) {
    //         let id = ident.unwrap();
    //         let id_str = id.to_string();
    //         let schema = avro_schema::from_syn(&id_str, ty);
    //         fids.push(id);
    //         ftys.push(ty.clone());
    //         fschemas.push(schema);
    //         fpositions.push(pos);
    //     }
    // }
    //
    // quote!(
    //     impl#impl_generics crate::utils::avro::Schematized for #id #ty_generics #where_clause {
    //         // Hello world! #impl_generics
    //         fn schematize(parent_space: Option<String>) -> avro_rs::schema::Schema {
    //             let new_namespace = if let Some(pns) = parent_space.as_ref() {
    //                 Some(vec![String::from(pns), String::from(#id_string)].join("."))
    //             } else {
    //                 Some(String::from(#id_string))
    //             };
    //             avro_rs::schema::Schema::Array()
    //         }
    //     }
    // )
}

#[proc_macro_derive(Schematize)]
pub fn derive(
    input: proc_macro::TokenStream,
) -> proc_macro::TokenStream {
    let DeriveInput {
        ident: id,
        data,
        generics,
        ..
    } = parse_macro_input!(input as DeriveInput);
    let r = match data {
        Data::Struct(DataStruct {
            fields: Fields::Named(nfs),
            ..
        }) => gen_struct_impl(id, generics, nfs),
        Data::Enum(DataEnum { variants, .. }) => gen_enum_impl(
            id,
            generics,
            variants
                .iter()
                .map(|v| v.clone())
                .collect::<Vec<Variant>>(),
        ),
        _ => panic!(
            "Can only handle structs with named fields, and enums."
        ),
    };

    proc_macro::TokenStream::from(r)
}
