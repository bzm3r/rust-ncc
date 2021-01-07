use crate::world::{Cells, Snapshot};
use bincode::serialize_into;
use bson::Document;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::path::PathBuf;

#[derive(Clone, Copy)]
pub enum Format {
    Bson,
    Json,
    Bincode,
}

pub fn get_file(
    out_dir: &PathBuf,
    title: &str,
    format: Format,
) -> io::Result<File> {
    let mut path = out_dir.clone();

    let ext = match format {
        Format::Bson => ".bson".to_string(),
        Format::Json => ".json".to_string(),
        Format::Bincode => ".binc".to_string(),
    };

    #[cfg(features = "validation")]
    path.push(format!("history_dbg_{}.{}", title, ext));
    #[cfg(not(features = "validation"))]
    path.push(format!("history_{}.{}", title, ext));

    OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&path)
}

pub fn save_compact(
    data: Vec<(u32, Cells)>,
    out_dir: &PathBuf,
    formats: Vec<Format>,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    for format in formats {
        let mut f = get_file(out_dir, title, format)?;

        match format {
            Format::Bson => {
                let mut doc = Document::new();
                doc.insert("data", bson::to_bson(&data)?);
                doc.to_writer(&mut f);
            }
            Format::Json => serde_json::to_writer(&mut f, &data)?,
            Format::Bincode => serialize_into(&mut f, &data)?,
        };
    }

    Ok(())
}

pub fn save_full(
    data: Vec<(u32, Snapshot)>,
    out_dir: &PathBuf,
    formats: Vec<Format>,
    title: &str,
) -> Result<(), Box<dyn Error>> {
    for format in formats {
        let mut f = get_file(out_dir, title, format)?;

        match format {
            Format::Bson => {
                let mut doc = Document::new();
                let bson = bson::to_bson(&data)?;
                doc.insert("data", bson).unwrap();
                doc.to_writer(&mut f);
            }
            Format::Json => serde_json::to_writer(&mut f, &data)?,
            Format::Bincode => serialize_into(&mut f, &data)?,
        };
    }

    Ok(())
}
