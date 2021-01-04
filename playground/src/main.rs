mod utils;

use avro_rs::schema::{
    Name, RecordField, RecordFieldOrder, UnionSchema,
};
use avro_rs::{to_value, Schema, Writer};
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Copy, Serialize, Deserialize, Schematize)]
struct MyStruct {
    x: bool,
}

#[derive(Serialize, Deserialize)]
enum MyEnum {
    Zero(bool),
    One(MyStruct),
    Two(MyStruct, bool),
    None,
}

fn main() {
    // let struct_schema = Schema::Record {
    //     name: Name {
    //         name: "MyStruct".to_string(),
    //         namespace: None,
    //         aliases: None,
    //     },
    //     doc: None,
    //     fields: vec![RecordField {
    //         name: "x".to_string(),
    //         doc: None,
    //         default: None,
    //         schema: Schema::Boolean,
    //         order: RecordFieldOrder::Ascending,
    //         position: 0,
    //     }],
    //     lookup: {
    //         let mut m: HashMap<String, usize> = HashMap::new();
    //         m.insert("x".to_string(), 0);
    //         m
    //     },
    // };
    //
    // let zero_schema = Schema::Record {
    //     name: Name {
    //         name: "Zero".to_string(),
    //         namespace: Some("MyEnum".to_string()),
    //         aliases: None,
    //     },
    //     doc: None,
    //     fields: vec![RecordField {
    //         name: "0".to_string(),
    //         doc: None,
    //         default: None,
    //         schema: Schema::Boolean,
    //         order: RecordFieldOrder::Ascending,
    //         position: 0,
    //     }],
    //     lookup: {
    //         let mut m: HashMap<String, usize> = HashMap::new();
    //         m.insert("0".to_string(), 0);
    //         m
    //     },
    // };
    // let one_schema = Schema::Record {
    //     name: Name {
    //         name: "One".to_string(),
    //         namespace: Some("MyEnum".to_string()),
    //         aliases: None,
    //     },
    //     doc: None,
    //     fields: vec![RecordField {
    //         name: "0".to_string(),
    //         doc: None,
    //         default: None,
    //         schema: struct_schema.clone(),
    //         order: RecordFieldOrder::Ascending,
    //         position: 0,
    //     }],
    //     lookup: {
    //         let mut m: HashMap<String, usize> = HashMap::new();
    //         m.insert("0".to_string(), 0);
    //         m
    //     },
    // };
    // let two_schema = Schema::Record {
    //     name: Name {
    //         name: "Two".to_string(),
    //         namespace: Some("MyEnum".to_string()),
    //         aliases: None,
    //     },
    //     doc: None,
    //     fields: vec![
    //         RecordField {
    //             name: "0".to_string(),
    //             doc: None,
    //             default: None,
    //             schema: struct_schema.clone(),
    //             order: RecordFieldOrder::Ascending,
    //             position: 0,
    //         },
    //         RecordField {
    //             name: "1".to_string(),
    //             doc: None,
    //             default: None,
    //             schema: Schema::Boolean,
    //             order: RecordFieldOrder::Ascending,
    //             position: 1,
    //         },
    //     ],
    //     lookup: {
    //         let mut m: HashMap<String, usize> = HashMap::new();
    //         m.insert("0".to_string(), 0);
    //         m
    //     },
    // };
    // let none_schema = Schema::Null;
    //
    // let enum_schema = Schema::Union(
    //     UnionSchema::new(vec![
    //         none_schema,
    //         zero_schema,
    //         one_schema,
    //         two_schema,
    //     ])
    //     .unwrap(),
    // );

    // let mut struct_writer = Writer::new(&struct_schema, Vec::new());
    // the structure models our Record schema
    let ms = MyStruct { x: false };
    let me = MyEnum::Two(ms.clone(), true);
    //
    // // schema validation happens here
    // struct_writer.append_ser(ms).unwrap();
    //
    // // this is how to get back the resulting avro bytecode
    // // this performs a flush operation to make sure data is written, so it can fail
    // // you can also call `writer.flush()` yourself without consuming the writer
    // let encoded = struct_writer.into_inner();

    let v = to_value(me).unwrap();
    println!("{:?}", v);
}
