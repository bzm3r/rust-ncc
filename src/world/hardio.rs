use crate::world::{DeepHistory, History};
use bincode::{deserialize_from, serialize_into};
use std::error::Error;
use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::path::PathBuf;

#[allow(unused)]
#[derive(Clone, Copy)]
pub enum Format {
    Cbor,
    Json,
    Bincode,
}

pub fn get_file_name(
    compact: bool,
    format: Format,
    title: &str,
) -> String {
    let ext = match format {
        Format::Cbor => "cbor".to_string(),
        Format::Json => "json".to_string(),
        Format::Bincode => "binc".to_string(),
    };

    let dat_ty = if compact {
        "compact".to_string()
    } else {
        "full".to_string()
    };

    if cfg!(features = "validate") {
        format!("history_validated_{}_{}.{}", dat_ty, title, ext)
    } else {
        format!("history_{}_{}.{}", dat_ty, title, ext)
    }
}

pub fn get_file_write(
    compact: bool,
    out_dir: &PathBuf,
    title: &str,
    format: Format,
) -> io::Result<File> {
    let mut path = out_dir.clone();
    path.push(get_file_name(compact, format, title));

    OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&path)
}

pub fn read_file(
    compact: bool,
    out_dir: &PathBuf,
    title: &str,
    format: Format,
) -> io::Result<File> {
    let mut path = out_dir.clone();
    path.push(get_file_name(compact, format, title));

    OpenOptions::new().read(true).open(&path)
}

pub fn save_compact(
    data: History,
    out_dir: &PathBuf,
    formats: Vec<Format>,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    for format in formats {
        let mut f = get_file_write(true, out_dir, title, format)?;

        match format {
            Format::Cbor => {
                serde_cbor::to_writer(&mut f, &data)?;
            }
            Format::Json => serde_json::to_writer(&mut f, &data)?,
            Format::Bincode => serialize_into(&mut f, &data)?,
        };
    }

    Ok(())
}

pub fn save_full(
    data: DeepHistory,
    out_dir: &PathBuf,
    formats: Vec<Format>,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    for format in formats {
        let mut f = get_file_write(false, out_dir, title, format)?;

        match format {
            Format::Cbor => {
                serde_cbor::to_writer(&mut f, &data)?;
            }
            Format::Json => serde_json::to_writer(&mut f, &data)?,
            Format::Bincode => serialize_into(&mut f, &data)?,
        };
    }

    Ok(())
}

pub fn load_compact(
    out_dir: &PathBuf,
    format: Format,
    title: &str,
) -> History {
    let mut f = read_file(true, out_dir, title, format).unwrap();

    match format {
        Format::Cbor => serde_cbor::from_reader(&mut f).unwrap(),
        Format::Json => serde_json::from_reader(&mut f).unwrap(),
        Format::Bincode => deserialize_from(&mut f).unwrap(),
    }
}

pub fn load_binc_from_path(file_path: &PathBuf) -> History {
    if let Some(ext) = file_path.extension() {
        match ext.to_str().unwrap() {
            "binc" => {
                let mut f = OpenOptions::new()
                    .read(true)
                    .open(&file_path)
                    .unwrap();
                deserialize_from(&mut f).unwrap()
            }
            _ => panic!(
                "file path has unknown extension: {}",
                ext.to_str().unwrap()
            ),
        }
    } else {
        panic!(
            "no file extension in path: {}",
            file_path.to_str().unwrap()
        )
    }
}

pub fn load_full(
    out_dir: &PathBuf,
    format: Format,
    title: &str,
) -> DeepHistory {
    let mut f = read_file(false, out_dir, title, format).unwrap();

    match format {
        Format::Cbor => serde_cbor::from_reader(&mut f).unwrap(),
        Format::Json => serde_json::from_reader(&mut f).unwrap(),
        Format::Bincode => deserialize_from(&mut f).unwrap(),
    }
}
