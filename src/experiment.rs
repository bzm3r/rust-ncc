// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use serde::Deserialize;
use crate::parameters::UserParamsOverrides;
use crate::math::P2D;
use std::path::Path;
use crate::quantity::Length;

#[derive(Deserialize)]
pub struct GroupLayout {
    pub width: u32,
    pub height: u32,
    pub centroid: [Length; 2],
}

#[derive(Deserialize)]
pub struct CellGroup {
    pub num_cells: u32,
    pub layout: GroupLayout,
    pub parameter_overrides: UserParamsOverrides,
}

#[derive(Deserialize)]
pub struct Experiment {
    pub title: String,
    pub date: String,
    pub cell_groups: Vec<CellGroup>,
}

pub fn load_experiment(name: &str) -> Experiment {
    let path_string = format!("./experiments/{}.toml", name);
    let path = Path::new(&path_string);

    let err_string = format!("could not read from file: {}", path_string);
    let out = std::fs::read_to_string(path).expect(&err_string);

    toml::from_str(&out).unwrap()
}
