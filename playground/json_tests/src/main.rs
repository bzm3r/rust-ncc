use serde::{Deserialize, Serialize};
use std::fs::OpenOptions;
use std::io::Read;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum DistribDef {
    Random { frac: f64 },
    SpecificRandom { frac: f64, marked_verts: Vec<usize> },
    SpecificUniform { frac: f64, marked_verts: Vec<usize> },
}

impl Default for DistribDef {
    fn default() -> Self {
        DistribDef::Random { frac: 0.1 }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RgtpDistribDef {
    acts: DistribDef,
    inacts: DistribDef,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RgtpDistribDefs {
    rac: RgtpDistribDef,
    rho: RgtpDistribDef,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct PairRgtpDistribDefs {
    cell0: RgtpDistribDefs,
    cell1: RgtpDistribDefs,
}

#[allow(clippy::large_enum_variant)]
#[derive(Clone, Deserialize, Serialize, Debug)]
pub enum ExperimentType {
    PyCompare {
        num_cells: usize,
    },
    Pair {
        sep_in_cell_diams: usize,
        rgtp_distrib_defs_per_cell: PairRgtpDistribDefs,
    },
    NCells {
        num_cells: usize,
    },
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub struct EulerOpts {
    pub num_int_steps: usize,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub struct Rkdp5Opts {
    pub max_iters: usize,
    pub atol: f64,
    pub rtol: f64,
    pub init_h_scale: f64,
}

impl Default for Rkdp5Opts {
    fn default() -> Self {
        Rkdp5Opts {
            max_iters: 20,
            atol: 1e-3,
            rtol: 1e-3,
            init_h_scale: 0.1,
        }
    }
}

#[derive(Clone, Copy, Deserialize, Serialize, PartialEq, Debug)]
pub enum IntegratorOpts {
    Euler(EulerOpts),
    EulerDebug(EulerOpts),
    Rkdp5(Rkdp5Opts),
}

#[derive(Clone, Copy, Deserialize, Serialize, PartialEq, Debug)]
pub struct AnimationOptions {
    label_verts: bool,
    label_cells: bool,
    follow_group: bool,
    show_trails: bool,
    rgtp_scale: f64,
}

#[derive(Clone, Deserialize, Serialize, Debug)]
pub struct ParsedExpArgs {
    description: String,
    experiment_type: ExperimentType,
    final_t: f64,
    cil: f64,
    coa: Option<f64>,
    adh: Option<f64>,
    cal: Option<f64>,
    snap_period: f64,
    randomization: bool,
    rgtp_distrib_defs: Option<RgtpDistribDefs>,
    seeds: Vec<u64>,
    int_opts: IntegratorOpts,
    ani_opts: Vec<AnimationOptions>,
}

fn main() {
    let mut test_json_file = OpenOptions::new()
        .read(true)
        .open("./pair_sep.json")
        .unwrap();
    let mut test_json_out = String::new();
    test_json_file.read_to_string(&mut test_json_out).unwrap();
    let parsed_args: ParsedExpArgs =
        serde_json::from_str(&test_json_out).unwrap();
    println!("{:#?}", parsed_args);
}
