use crate::exp_setup::{ExperimentType, RgtpDistribDefs};
use crate::parameters::quantity::Time;
use crate::world::{EulerOpts, IntegratorOpts, RkOpts};
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::error;
use std::fs::OpenOptions;
use std::io::Read;
use std::path::PathBuf;

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
struct ParsedRkOpts {
    max_iters: Option<usize>,
    atol: Option<f64>,
    rtol: Option<f64>,
    init_h_scale: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
struct ParsedEulerOpts {
    num_int_steps: Option<usize>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
enum ParsedIntOpts {
    RkDp5(ParsedRkOpts),
    Euler(ParsedEulerOpts),
    EulerDebug(ParsedEulerOpts),
}

impl Default for ParsedIntOpts {
    fn default() -> Self {
        ParsedIntOpts::RkDp5(ParsedRkOpts::default())
    }
}

impl From<ParsedIntOpts> for IntegratorOpts {
    fn from(parsed_opts: ParsedIntOpts) -> Self {
        match parsed_opts {
            ParsedIntOpts::RkDp5(opts) => {
                let default = RkOpts::default();
                IntegratorOpts::Rkdp5(RkOpts {
                    max_iters: opts
                        .max_iters
                        .unwrap_or(default.max_iters),
                    atol: opts.atol.unwrap_or(default.atol),
                    rtol: opts.rtol.unwrap_or(default.rtol),
                    init_h_scale: opts
                        .init_h_scale
                        .unwrap_or(default.init_h_scale),
                })
            }
            ParsedIntOpts::Euler(opts) => {
                IntegratorOpts::Euler(EulerOpts {
                    num_int_steps: opts.num_int_steps.unwrap_or(10),
                })
            }
            ParsedIntOpts::EulerDebug(opts) => {
                IntegratorOpts::EulerDebug(EulerOpts {
                    num_int_steps: opts.num_int_steps.unwrap_or(10),
                })
            }
        }
    }
}

#[derive(Clone, Copy, Deserialize, Serialize, Debug)]
pub struct AnimationOptions {
    label_verts: bool,
    label_cells: bool,
    follow_group: bool,
    show_trails: bool,
    arrow_scale: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
struct ParsedExpArgs {
    ty: ExperimentType,
    final_t: f64,
    cil_mag: f64,
    coa_mag: Option<f64>,
    adh_mag: Option<f64>,
    cal_mag: Option<f64>,
    snap_period: f64,
    max_on_ram: Option<usize>,
    randomization: bool,
    rgtp_distrib_defs: Option<RgtpDistribDefs>,
    seeds: Vec<u64>,
    int_opts: ParsedIntOpts,
    ani_opts: Vec<AnimationOptions>,
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct ExperimentArgs {
    pub file_name: String,
    pub ty: ExperimentType,
    pub final_t: Time,
    pub cil_mag: f64,
    pub coa_mag: Option<f64>,
    pub cal_mag: Option<f64>,
    pub adh_mag: Option<f64>,
    pub snap_period: Time,
    pub max_on_ram: usize,
    pub randomization: bool,
    pub seeds: Vec<u64>,
    pub int_opts: IntegratorOpts,
    pub rgtp_distrib_defs: RgtpDistribDefs,
}

impl TryFrom<&PathBuf> for ExperimentArgs {
    type Error = Box<dyn error::Error>;

    fn try_from(
        json_path: &PathBuf,
    ) -> Result<ExperimentArgs, Box<dyn error::Error>> {
        let mut f = OpenOptions::new().read(true).open(json_path)?;
        let mut json_out = String::new();
        f.read_to_string(&mut json_out)?;
        let ParsedExpArgs {
            ty,
            final_t,
            cil_mag,
            coa_mag,
            adh_mag,
            cal_mag,
            snap_period,
            max_on_ram,
            randomization,
            rgtp_distrib_defs,
            seeds,
            int_opts,
            ..
        } = serde_json::from_str(&json_out).unwrap();
        let file_name: String = json_path
            .file_stem()
            .unwrap_or_else(|| {
                panic!(
                    "could not extract file name from path: {}",
                    &json_path.to_str().unwrap()
                )
            })
            .to_str()
            .unwrap_or_else(|| {
                panic!(
                    "could not convert file name to string from path: {}",
                    &json_path.to_str().unwrap()
                )
            })
            .into();
        let exp_args = ExperimentArgs {
            file_name,
            ty: ty.clone(),
            final_t: Time(final_t),
            cil_mag,
            coa_mag,
            cal_mag,
            adh_mag,
            snap_period: Time(snap_period),
            max_on_ram: max_on_ram.unwrap_or(1000),
            randomization,
            seeds,
            int_opts: int_opts.into(),
            rgtp_distrib_defs: rgtp_distrib_defs
                .unwrap_or(RgtpDistribDefs::default()),
        };
        Ok(exp_args)
    }
}
