#![allow(clippy::too_many_arguments)]
#![allow(clippy::large_enum_variant)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::suspicious_arithmetic_impl)]

pub mod cell;
pub mod exp_setup;
pub mod hardio;
pub mod interactions;
pub mod math;
pub mod parameters;
pub mod utils;
pub mod world;

use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::error;
use std::fs::{create_dir_all, OpenOptions};
use std::io::Read;
use std::path::PathBuf;

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
    type Error = Box<dyn error::Error>;

    fn try_from(
        json_path: &PathBuf,
    ) -> Result<Self, Box<dyn error::Error>> {
        let mut f = OpenOptions::new().read(true).open(json_path)?;
        let mut json_out = String::new();
        f.read_to_string(&mut json_out)?;
        let directories = serde_json::from_str(&json_out)?;
        Ok(directories)
    }
}
