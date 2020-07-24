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
    b: SubState,
}

impl State {
    fn random() -> State {
        State {
            a: P2D::random_vec(4),
            b: SubState::random(),
        }
    }
}

fn save_schema(name: &str, raw_schema: &str, output_dir: &PathBuf) {
    let mut avsc_path = output_dir.clone();
    avsc_path.push(format!("{}_schema", name));
    avsc_path.set_extension("avsc");

    let mut f = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(avsc_path)
        .unwrap();
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
    let raw_schema = State::raw_schema();
    save_schema(name, &raw_schema, output_dir);
    let schema = Schema::parse_str(&raw_schema).unwrap();
    let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
    writer.append_ser(state).unwrap();
    let encoded = writer.into_inner().unwrap();
    save_data(name, &encoded, output_dir);
}

fn main() {
    let state = State::random();
    let output_dir = PathBuf::from("C:\\Users\\bhmer\\Desktop\\rust-ncc\\playground\\output");
    save_state(&state, &output_dir);
    println!("{:?}", state);
}
