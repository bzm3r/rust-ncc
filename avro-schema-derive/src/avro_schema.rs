use quote::quote;

pub(crate) fn from_syn(fstr: &str, ty: &syn::Type) -> proc_macro2::TokenStream {
    match ty {
        syn::Type::Path(tp) => from_path(fstr, tp),
        syn::Type::Array(ta) => {
            let inner = from_syn(fstr, &ta.elem);
            quote!(avro_rs::schema::Schema::Array(Box::new(#inner)))
        }
        _ => panic!("Schematize: cannot handle non-Path or Array syn::Type."),
    }
}

pub fn process_id(fstr: &str, id: &syn::Ident) -> proc_macro2::TokenStream {
    let id_string = id.to_string();
    match id_string.as_str() {
        "bool" => quote!(avro_rs::schema::Schema::Boolean),
        "i32" | "u32" => quote!(avro_rs::schema::Schema::Int),
        "i64" | "u64" => quote!(avro_rs::schema::Schema::Long),
        "f32" => quote!(avro_rs::schema::Schema::Float),
        "f64" => quote!(avro_rs::schema::Schema::Double),
        "String" => quote!(avro_rs::schema::Schema::String),
        _ => {
            quote!(#id::schematize(Some(vec![new_namespace.clone().unwrap(), String::from(#fstr)].join("."))))
        }
    }
}

fn process_single_seg_path(fstr: &str, seg: &syn::PathSegment) -> proc_macro2::TokenStream {
    let seg_id_string = seg.ident.to_string();
    match seg_id_string.as_str() {
        "Vec" => match &seg.arguments {
            syn::PathArguments::AngleBracketed(angle_args) => {
                let args = &angle_args.args;
                match args.len() {
                        1 => match args.first().unwrap() {
                            syn::GenericArgument::Type(t) => {
                                let inner = from_syn(fstr, t);
                                quote!(avro_rs::schema::Schema::Array(Box::new(#inner)))
                            },
                            _ => panic!("Schematize: encountered variant of syn::GenericArgument other than Type."),
                        },
                        _ => panic!("Schematize: cannot handle multi-arg vecs."),
                    }
            }
            _ => panic!("Schematize: encountered Vec without <>."),
        },
        "Box" => match &seg.arguments {
            syn::PathArguments::AngleBracketed(angle_args) => {
                let args = &angle_args.args;
                match args.len() {
                        1 => match args.first().unwrap() {
                            syn::GenericArgument::Type(t) => from_syn(fstr, t),
                            _ => panic!("Schematize: encountered variant of syn::GenericArgument other than Type."),
                        },
                        _ => panic!("Schematize: cannot handle multi-arg boxes."),
                    }
            }
            _ => panic!("Schematize: encountered Box without <>."),
        },
        _ => quote!("Schematize: encountered single-seg path that is not Box or Vec."),
    }
}

fn from_path(fstr: &str, tp: &syn::TypePath) -> proc_macro2::TokenStream {
    if tp.path.segments.is_empty() {
        panic!("Schematize: path contains no segments.");
    } else if let Some(id) = tp.path.get_ident() {
        process_id(fstr, id)
    } else if tp.path.segments.len() == 1 {
        process_single_seg_path(fstr, tp.path.segments.first().unwrap())
    } else {
        let ids = tp
            .path
            .segments
            .iter()
            .filter_map(|seg| match seg.arguments {
                syn::PathArguments::None => Some(seg.ident.clone()),
                _ => None,
            })
            .collect::<Vec<syn::Ident>>();
        quote!(#(#ids)::*)
    }
}
