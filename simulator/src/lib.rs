pub mod cell;
pub mod exp_setup;
pub mod hardio;
pub mod interactions;
pub mod math;
pub mod parameters;
pub mod toml_parse;
pub mod utils;
pub mod world;

use crate::toml_parse::{
    get_value, parse_path, parse_table, ParseErr,
};
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::fs::{create_dir_all, OpenOptions};
use std::io::Read;
use std::path::PathBuf;
use toml::Value;

/// Number of vertices per model cell.
pub const NVERTS: usize = 16;
/// Default directory where simulation output will be placed.
pub const DEFAULT_OUTPUT_DIR: &str = "B:\\rust-ncc\\output";

#[derive(Clone, Deserialize, Serialize)]
pub struct Directories {
    pub out: PathBuf,
    pub exp: PathBuf,
    py_main: PathBuf,
}

impl Directories {
    pub fn make(&self) {
        [&self.out, &self.exp]
            .iter()
            .for_each(|p| create_dir_all(p).unwrap());
    }
}

impl TryFrom<&PathBuf> for Directories {
    type Error = ParseErr;

    fn try_from(toml_path: &PathBuf) -> Result<Self, Self::Error> {
        let mut f = OpenOptions::new()
            .read(true)
            .open(toml_path)
            .map_err(|e| {
                ParseErr::FileOpen(format!(
                    "{:?}: {}",
                    e,
                    toml_path.to_str().unwrap()
                ))
            })?;
        let raw_cfg = {
            let mut out = String::new();
            let _ = f.read_to_string(&mut out).map_err(|_| {
                ParseErr::FileParse(String::from(
                    toml_path.to_str().unwrap(),
                ))
            })?;
            let value = out.parse::<Value>().map_err(|_| {
                ParseErr::FileParse(String::from(
                    toml_path.to_str().unwrap(),
                ))
            })?;
            parse_table(&value)?
        };
        let out = {
            let value = get_value(&raw_cfg, "out_dir")?;
            parse_path(&value)?
        };
        let exp = {
            let value = get_value(&raw_cfg, "exp_dir")?;
            parse_path(&value)?
        };
        let py_main = {
            let value = get_value(&raw_cfg, "py_main")?;
            parse_path(&value)?
        };
        Ok(Directories { out, exp, py_main })
    }
}
