use crate::world::Cells;
use avro_rs::{from_value, Reader, Schema};
use std::fs::{read, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn save_schema(name: &str, schema: &Schema, output_dir: &PathBuf) {
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

pub fn load_history(history_path: &Path) -> Vec<Cells> {
    let mut r = vec![];
    let history = read(&history_path).unwrap();
    let reader = Reader::new(&history[..]).unwrap();

    // value is a Result in case the read operation fails
    for value in reader {
        r.push(from_value::<Cells>(&value.unwrap()).unwrap())
    }

    r
}
