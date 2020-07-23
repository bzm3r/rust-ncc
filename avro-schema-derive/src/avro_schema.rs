use quote::{format_ident, quote};

use crate::print;
use std::borrow::ToOwned;
use std::boxed::Box;

pub enum AvroType {
    Primitive(String),
    Record(String),
    Array(Box<AvroType>),
}

pub enum RustType {
    Primitive(String),
    Path(PathType),
    Array(ArrayType),
    Vector(VectorType),
    Box(BoxType),
}

impl RustType {
    pub fn from_syn_type(ty: &syn::Type) -> RustType {
        match ty {
            syn::Type::Path(tp) => PathType::process(&tp),
            syn::Type::Array(ta) => RustType::Array(ArrayType::process(ta)),
            _ => panic!("Avro schema generator: cannot handle non-Path or Array syn::Type."),
        }
    }

    pub fn try_primitive(id: &syn::Ident) -> Option<RustType> {
        let id_string = id.to_string();
        match id_string.as_str() {
            "bool" | "i32" | "u32" | "i64" | "f32" | "f64" | "String" => {
                Some(RustType::Primitive(id_string))
            }
            _ => None,
        }
    }
}

pub struct ArrayType {
    inner: Box<RustType>,
}

impl ArrayType {
    fn process(ta: &syn::TypeArray) -> ArrayType {
        ArrayType {
            inner: Box::new(RustType::from_syn_type(&ta.elem)),
        }
    }
}

pub struct VectorType {
    inner: Box<RustType>,
}

impl VectorType {
    fn process(angle_args: &syn::AngleBracketedGenericArguments) -> RustType {
        let args = &angle_args.args;

        if args.len() == 1 {
            let inner = match args.first().unwrap() {
                syn::GenericArgument::Type(t) => Box::new(RustType::from_syn_type(t)),
                _ => panic!("VecType: do not know how to handle non-Type GenericArguments."),
            };
            RustType::Vector(VectorType { inner })
        } else {
            panic!(
                "Don't know how to handle Vec{}",
                print::format_angled_args(angle_args)
            )
        }
    }
}

pub struct BoxType {
    inner: Box<RustType>,
}

impl BoxType {
    fn process(angle_args: &syn::AngleBracketedGenericArguments) -> RustType {
        let args = &angle_args.args;

        if args.len() == 1 {
            let inner = match args.first().unwrap() {
                syn::GenericArgument::Type(t) => Box::new(RustType::from_syn_type(t)),
                _ => panic!("BoxType: do not know how to handle non-Type GenericArguments."),
            };
            RustType::Box(BoxType { inner })
        } else {
            panic!(
                "Don't know how to handle Box{}",
                print::format_angled_args(angle_args)
            )
        }
    }
}

pub struct PathType {
    idents: Vec<String>,
}

impl PathType {
    fn process(tp: &syn::TypePath) -> RustType {
        if tp.path.segments.is_empty() {
            panic!("Avro schema: cannot not handle TypePath with Path containing no segments.");
        } else if let Some(id) = tp.path.get_ident() {
            if let Some(primitive) = RustType::try_primitive(id) {
                primitive
            } else {
                RustType::Path(PathType {
                    idents: vec![id.to_string()],
                })
            }
        } else if tp.path.segments.len() == 1 {
            let first = tp.path.segments.first().unwrap();
            let id = first.ident.to_string();
            if &id == "Vec" {
                match &first.arguments {
                    syn::PathArguments::AngleBracketed(angle_args) => {
                        VectorType::process(angle_args)
                    }
                    syn::PathArguments::None => panic!(
                        "Encountered a Vec with no angle arguments: {}",
                        print::format_path(tp)
                    ),
                    syn::PathArguments::Parenthesized(_) => panic!(
                        "Encountered a Vec with paranthesis arguments: {}",
                        print::format_path(tp)
                    ),
                }
            } else if &id == "Box" {
                match &first.arguments {
                    syn::PathArguments::AngleBracketed(angle_args) => BoxType::process(angle_args),
                    syn::PathArguments::None => panic!(
                        "Encountered a Box with no angle arguments: {}",
                        print::format_path(tp)
                    ),
                    syn::PathArguments::Parenthesized(_) => panic!(
                        "Encountered a Box with paranthesis arguments: {}",
                        print::format_path(tp)
                    ),
                }
            } else {
                panic!(
                    "Do not know how to handle TypePath: {}",
                    print::format_path(tp)
                )
            }
        } else {
            RustType::Path(PathType {
                idents: tp
                    .path
                    .segments
                    .iter()
                    .map(|seg| match PathType::process_seg(seg) {
                        Some(id) => id,
                        None => panic!(
                            "Do not know how to handle TypePath: {}",
                            print::format_path(tp)
                        ),
                    })
                    .collect(),
            })
        }
    }

    fn process_seg(seg: &syn::PathSegment) -> Option<String> {
        match seg.arguments {
            syn::PathArguments::None => Some(seg.ident.to_string()),
            _ => None,
        }
    }
}

impl AvroType {
    pub fn from_rust_type(ty: &RustType) -> AvroType {
        match ty {
            RustType::Vector(vt) => AvroType::Array(Box::new(AvroType::from_rust_type(&vt.inner))),
            RustType::Box(bt) => AvroType::from_rust_type(&bt.inner),
            RustType::Array(at) => AvroType::Array(Box::new(AvroType::from_rust_type(&at.inner))),
            RustType::Path(pt) => AvroType::Record(pt.idents.join("::")),
            RustType::Primitive(pt) => AvroType::Primitive(Self::primitive(pt).unwrap()),
        }
    }

    pub fn primitive(rust_primitive: &str) -> Option<String> {
        match rust_primitive {
            "bool" => Some("boolean".to_owned()),
            // TODO: a usize::MAX > isize::MAX > u32::MAX > i32::MAX...need to handle gracefully!
            "i32" | "u32" => Some("int".to_owned()),
            "i64" | "u64" => Some("long".to_owned()),
            "f32" => Some("float".to_owned()),
            "f64" => Some("double".to_owned()),
            "String" => Some("string".to_owned()),
            _ => None,
        }
    }

    pub fn schema_macro_entry(&self) -> proc_macro2::TokenStream {
        match self {
            AvroType::Record(r) => {
                let rid = format_ident!("{}", r);
                quote!(#rid::raw_schema())
            }
            AvroType::Primitive(p) => {
                let p = p.to_string();
                quote!(format!("\"{}\"", #p))
            }
            AvroType::Array(inner) => {
                let se = inner.schema_macro_entry();
                quote!(format!("{{ \"type\": \"array\", \"items\": {} }}", #se))
            }
        }
    }
}
