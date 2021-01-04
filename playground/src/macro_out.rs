#![feature(prelude_import)]
#[prelude_import]
use std::prelude::v1::*;
#[macro_use]
extern crate std;
use avro_rs::Schema;
use avro_schema_derive::Schematize;

pub trait Schematize {
    fn schematize(
        parent_space: Option<String>,
    ) -> avro_rs::schema::Schema;
}
impl Schematize for bool {
    fn schematize(_parent_space: Option<String>) -> Schema {
        avro_rs::schema::Schema::Boolean
    }
}
struct MyStruct<T: Schematize> {
    a: T,
}
impl<T: Schematize> Schematize for MyStruct<T> {
    fn schematize(
        parent_space: Option<String>,
    ) -> avro_rs::schema::Schema {
        let new_namespace = if let Some(pns) = parent_space.as_ref() {
            Some(
                <[_]>::into_vec(box [
                    String::from(pns),
                    String::from("MyStruct"),
                ])
                .join("."),
            )
        } else {
            Some(String::from("MyStruct"))
        };
        avro_rs::schema::Schema::Record {
            name: avro_rs::schema::Name {
                name: String::from("MyStruct"),
                namespace: parent_space.clone(),
                aliases: None,
            },
            doc: None,
            fields: <[_]>::into_vec(box [
                avro_rs::schema::RecordField {
                    name: std::string::String::from("a"),
                    doc: None,
                    default: None,
                    schema: T::schematize(
                        Some(
                            <[_]>::into_vec(box [
                                new_namespace.clone().unwrap(),
                                String::from("a"),
                            ])
                            .join("."),
                        ),
                        None,
                    ),
                    order: avro_rs::schema::RecordFieldOrder::Ignore,
                    position: 0usize,
                },
            ]),
            lookup: {
                let mut r: std::collections::HashMap<String, usize> =
                    std::collections::HashMap::new();
                r.insert(std::string::String::from("a"), 0usize);
                r
            },
        }
    }
}
struct MyOtherStruct {
    b: MyStruct<bool>,
}
impl Schematize for MyOtherStruct {
    fn schematize(
        parent_space: Option<String>,
    ) -> avro_rs::schema::Schema {
        let new_namespace = if let Some(pns) = parent_space.as_ref() {
            Some(
                <[_]>::into_vec(box [
                    String::from(pns),
                    String::from("MyOtherStruct"),
                ])
                .join("."),
            )
        } else {
            Some(String::from("MyOtherStruct"))
        };
        avro_rs::schema::Schema::Record {
            name: avro_rs::schema::Name {
                name: String::from("MyOtherStruct"),
                namespace: parent_space.clone(),
                aliases: None,
            },
            doc: None,
            fields: <[_]>::into_vec(box [
                avro_rs::schema::RecordField {
                    name: std::string::String::from("b"),
                    doc: None,
                    default: None,
                    schema: { () },
                    order: avro_rs::schema::RecordFieldOrder::Ignore,
                    position: 0usize,
                },
            ]),
            lookup: {
                let mut r: std::collections::HashMap<String, usize> =
                    std::collections::HashMap::new();
                r.insert(std::string::String::from("b"), 0usize);
                r
            },
        }
    }
}
#[allow(dead_code)]
fn main() {
    let x = MyStruct { a: true };
    {
        ::std::io::_print(::core::fmt::Arguments::new_v1(
            &["Hello, world! ", "\n"],
            &match (&x.a,) {
                (arg0,) => [::core::fmt::ArgumentV1::new(
                    arg0,
                    ::core::fmt::Display::fmt,
                )],
            },
        ));
    };
}
#[main]
pub fn main() -> () {
    extern crate test;
    test::test_main_static(&[])
}
