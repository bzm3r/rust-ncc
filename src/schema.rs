use avro_rs::Schema;
use std::path::PathBuf;
use std::fs::OpenOptions;
use std::io::Write;

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

