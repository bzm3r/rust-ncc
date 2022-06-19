// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use std::path::Path;
use std::process::Command;
use tracing::info;

pub fn execute_py_model(
    out_dir: &Path,
    py_main: &Path,
    run_python: bool,
    name: &str,
    final_t: f64,
    snap_period: f64,
    num_cells: u32,
    cil: f64,
    coa: Option<f64>,
) {
    if run_python {
        let py_cmd: String = format!(
            "python {} --name {} --final_t {} --snap_period {} --num_cells {} --cil {} --coa {} \
        --out_dir {}",
            py_main.display(),
            name,
            final_t,
            snap_period,
            num_cells,
            cil,
            coa.unwrap_or(0.0),
            out_dir.display(),
        );
        info!("executing command: {}", py_cmd);
        if cfg!(target_os = "windows") {
            Command::new("cmd")
                .args(&["/C", &py_cmd])
                .status()
                .expect("failed to execute process")
        } else {
            Command::new("sh")
                .args(&["-c", &py_cmd])
                .status()
                .expect("failed to execute process")
        };
    }
}
