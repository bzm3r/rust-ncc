use std::borrow::ToOwned;
use syn::{AngleBracketedGenericArguments, ParenthesizedGenericArguments, Type};

pub fn format_angled_args(angle_args: &AngleBracketedGenericArguments) -> String {
    format!(
        "<{}>",
        angle_args
            .args
            .iter()
            .map(|ga| match ga {
                syn::GenericArgument::Type(t) => {
                    match t {
                        syn::Type::Path(tp2) => format_path(&tp2),
                        _ => "?".to_owned(),
                    }
                }
                _ => "_".to_owned(),
            })
            .collect::<Vec<String>>()
            .join(",")
    )
}

pub fn format_paren_args(paren_args: &ParenthesizedGenericArguments) -> String {
    format!(
        "({})",
        paren_args
            .inputs
            .iter()
            .map(|t| {
                match t {
                    syn::Type::Path(tp2) => format_path(&tp2),
                    _ => "?".to_owned(),
                }
            })
            .collect::<Vec<String>>()
            .join(",")
    )
}

pub fn format_path(tp: &syn::TypePath) -> String {
    let path_segs: Vec<String> = tp
        .path
        .segments
        .iter()
        .map(|seg| {
            format!(
                "{}{}",
                seg.ident.to_string(),
                match &seg.arguments {
                    syn::PathArguments::AngleBracketed(angles) => format_angled_args(angles),
                    syn::PathArguments::Parenthesized(parens) => format_paren_args(parens),
                    syn::PathArguments::None => "".to_owned(),
                }
            )
        })
        .collect::<Vec<String>>();

    path_segs.join("::")
}

#[allow(dead_code)]
pub fn format_syn_type(ty: &syn::Type) -> String {
    match ty {
        syn::Type::Array(_) => "Array".to_owned(),
        Type::BareFn(_) => "BareFn".to_owned(),
        Type::Group(_) => "Group".to_owned(),
        Type::ImplTrait(_) => "ImplTrait".to_owned(),
        Type::Infer(_) => "Infer".to_owned(),
        Type::Macro(_) => "Macro".to_owned(),
        Type::Never(_) => "Never".to_owned(),
        Type::Paren(_) => "Paren".to_owned(),
        Type::Path(path) => format_path(path),
        Type::Ptr(_) => "Ptr".to_owned(),
        Type::Reference(_) => "Reference".to_owned(),
        Type::Slice(_) => "Slice".to_owned(),
        Type::TraitObject(_) => "TraitObject".to_owned(),
        Type::Tuple(_) => "Tuple".to_owned(),
        Type::Verbatim(_) => "Verbatim".to_owned(),
        Type::__Nonexhaustive => "__Nonexhaustive".to_owned(),
    }
}
