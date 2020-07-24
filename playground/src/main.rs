use avro_schema_derive::Schematize;
use rand::distributions::Uniform;
use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};
use std::fs::OpenOptions;
use std::path::PathBuf;
use std::io::Write;
use avro_rs::{Schema, Writer};

#[derive(Debug, Default, Deserialize, Schematize, Serialize)]
pub struct P2D {
    x: f32,
    y: f32,
}

impl P2D {
    pub fn random_vec(size: usize) -> Vec<P2D> {
        let mut r = vec![];
        let distrib: Uniform<f32> = Uniform::new_inclusive(0.0, 1.0);
        let mut rng = thread_rng();

        for _ in 0..size {
            r.push(P2D {
                x: rng.sample(distrib),
                y: rng.sample(distrib),
            })
        }
        r
    }
}

#[derive(Debug, Default, Deserialize, Schematize, Serialize)]
pub struct SubState {
    c: Vec<P2D>,
}

impl SubState {
    fn random() -> SubState {
        SubState {
            c: P2D::random_vec(4),
        }
    }
}

#[derive(Debug, Default, Deserialize, Schematize, Serialize)]
pub struct State {
    a: Vec<P2D>,
    b: Vec<P2D>,
}

impl State {
    fn random() -> State {
        State {
            a: P2D::random_vec(4),
            b: P2D::random_vec(4),
        }
    }
}

fn save_schema(name: &str, schema: &Schema, output_dir: &PathBuf) {
    let mut avsc_path = output_dir.clone();
    avsc_path.push(format!("{}_schema", name));
    avsc_path.set_extension("avsc");

    let mut f = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(avsc_path)
        .unwrap();
    let raw_schema = serde_json::to_string_pretty(&schema).unwrap();
    println!("{}", &raw_schema);
    f.write_all(raw_schema.as_bytes()).unwrap();
}

pub fn save_data(name: &str, encoded: &[u8], output_dir: &PathBuf) {
    let mut path = output_dir.clone();
    path.push(format!("{}_dat", name));
    path.set_extension("avro");

    let mut f = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)
        .unwrap();
    f.write_all(&encoded).unwrap();
}

pub fn save_state(state: &State, output_dir: &PathBuf) {
    let name = "test";
    let schema = State::schematize(None);
    save_schema(name, &schema, &output_dir);
    let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
    writer.append_ser(state).unwrap();
    let encoded = writer.into_inner().unwrap();
    save_data(name, &encoded, output_dir);
}
//
// fn deduplicate_schema(schema: &serde_json::Value, mut known_schemas: Vec<String>) -> (Vec<String>, serde_json::Value) {
//     let schema_type = schema["type"].to_string().replace("\"", "");
//     if schema_type.as_str() == "record" {
//         let schema_name = schema["name"].to_string().replace("\"", "");
//         println!("beginning deduplication of schema {}, current known_schemas: {}", schema_name, format!("{:?}", &known_schemas));
//         if known_schemas.iter().any(|s| s == &schema_name) {
//             println!("found known record schema: {}, deduplicating...", schema_name);
//             let new_schema = serde_json::json!({
//                 "type": schema_name,
//             });
//             (known_schemas, new_schema)
//         } else {
//             println!("new record schema: {}", &schema_name);
//             known_schemas.push(schema_name.clone());
//             let mut dedup_fields: Vec<serde_json::Value> = vec![];
//             if let serde_json::Value::Array(orig_fields) = &schema["fields"] {
//                 for f in orig_fields {
//                     println!("analyzing field {:?}", f);
//                     let (kn, dedup_f) = deduplicate_schema(&f, known_schemas.clone());
//                     dedup_fields.push(dedup_f);
//                     known_schemas = kn;
//                 }
//                 let mut r = schema.clone();
//                 r["fields"] = serde_json::Value::Array(dedup_fields);
//                 (known_schemas, r)
//             } else {
//                 panic!("found non-vec fields!")
//             }
//         }
//     } else if schema_type.as_str() == "array" {
//         if let serde_json::Value::Array(items) = &schema["items"] {
//             let mut dedup_items: Vec<serde_json::Value> = vec![];
//             for i in items {
//                 println!("analyzing item of type: {}", i["type"]);
//                 let (kn, dedup_i) = deduplicate_schema(&i, known_schemas.clone());
//                 dedup_items.push(dedup_i);
//                 known_schemas = kn;
//             };
//             let mut r = schema.clone();
//             r["items"] = serde_json::Value::Array(dedup_items);
//             (known_schemas, r)
//         } else {
//             panic!("found non-vec array!")
//         }
//     } else {
//         let (known_schemas, r) = deduplicate_schema(&schema["type"], known_schemas);
//         (known_schemas, r)
//     }
// }

fn main() {
    let state = State::random();
    save_state(&state, &PathBuf::from("C:\\users\\bhmer\\desktop\\rust-ncc\\playground\\output"));
}
