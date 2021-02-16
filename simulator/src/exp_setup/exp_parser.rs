use crate::exp_setup::ExperimentType;
use crate::parameters::quantity::Time;
use crate::toml_parse::{
    get_optional_f64, get_optional_path, get_optional_usize,
    get_value, parse_array, parse_f64, parse_string, parse_table,
    parse_u64, ParseErr,
};
use crate::world::{EulerOpts, IntegratorOpts, Rkdp5Opts};
use std::convert::TryFrom;
use std::fs::OpenOptions;
use std::io::Read;
use std::path::PathBuf;
use toml::value::Table;
use toml::Value;

fn get_exp_type(
    raw_exp_args: &Table,
) -> Result<ExperimentType, ParseErr> {
    let exp_type_val = get_value(raw_exp_args, "experiment_type")?;
    let exp_type_str = parse_string(&exp_type_val)?;
    let exp_type = match exp_type_str.as_str() {
        "Pair" => ExperimentType::Pair {
            sep_in_cell_diams: {
                let v = get_value(raw_exp_args, "sep_in_cell_diams")?;
                (parse_u64(&v)?) as usize
            },
        },
        "PyCompare" => ExperimentType::PyCompare {
            num_cells: {
                let v = get_value(raw_exp_args, "num_cells")?;
                (parse_u64(&v)?) as usize
            },
            py_main: get_optional_path(raw_exp_args, "py_main")?,
        },
        _ => {
            return Err(ParseErr::UnknownExperiment(exp_type_str));
        }
    };
    Ok(exp_type)
}

fn parse_seed(value: &Value) -> Result<Option<u64>, ParseErr> {
    let inner = parse_table(value)?;
    match inner.get("seed") {
        Some(seed) => {
            let s = parse_u64(seed)?;
            Ok(Some(s))
        }
        None => Ok(None),
    }
}

fn get_seeds(
    raw_exp_args: &Table,
) -> Result<Vec<Option<u64>>, ParseErr> {
    let value = get_value(raw_exp_args, "seeds")?;
    let raw_seeds = parse_array(&value)?;
    let mut seeds: Vec<Option<u64>> =
        Vec::with_capacity(raw_seeds.len());
    let mut found_none = false;
    for raw_seed in raw_seeds {
        let seed = parse_seed(&raw_seed)?;
        if seed.is_none() && found_none {
            return Err(ParseErr::MultipleNoneSeeds);
        } else {
            found_none = seed.is_none();
        }
        seeds.push(seed);
    }
    Ok(seeds)
}

pub enum IntegratorType {
    Euler,
    EulerDebug,
    Rkdp5,
}

pub fn get_integrator_type(
    raw_int_opts: &Table,
) -> Result<IntegratorType, ParseErr> {
    let int_val = get_value(raw_int_opts, "integrator")
        .unwrap_or_else(|e| panic!("{:?}", e));
    let int_type = parse_string(&int_val)?;
    match int_type.as_str() {
        "Euler" => Ok(IntegratorType::Euler),
        "EulerDebug" => Ok(IntegratorType::EulerDebug),
        "Rkdp5" => Ok(IntegratorType::Rkdp5),
        _ => Err(ParseErr::UnknownIntegrator(int_type)),
    }
}

pub fn parse_euler_opts(
    raw_int_opts: &Table,
) -> Result<EulerOpts, ParseErr> {
    let default = EulerOpts::default();
    let num_int_steps =
        get_optional_usize(raw_int_opts, "num_int_steps")?
            .unwrap_or(default.num_int_steps);
    Ok(EulerOpts { num_int_steps })
}

pub fn parse_rkdp5_opts(
    raw_int_opts: &Table,
) -> Result<Rkdp5Opts, ParseErr> {
    let default = Rkdp5Opts::default();
    let max_iters = get_optional_usize(raw_int_opts, "max_iters")?
        .unwrap_or(default.max_iters);
    let atol = get_optional_f64(raw_int_opts, "atol")?
        .unwrap_or(default.atol);
    let rtol = get_optional_f64(raw_int_opts, "rtol")?
        .unwrap_or(default.rtol);
    let init_h_scale =
        get_optional_f64(raw_int_opts, "init_h_scale")?
            .unwrap_or(default.init_h_scale);
    Ok(Rkdp5Opts {
        max_iters,
        atol,
        rtol,
        init_h_scale,
    })
}

pub fn get_integrator_opts(raw_exp_args: &Table) -> IntegratorOpts {
    let value = get_value(raw_exp_args, "int_opts").unwrap();
    let raw_int_opts = parse_table(&value).unwrap();
    let int_type = get_integrator_type(&raw_int_opts).unwrap();
    match int_type {
        IntegratorType::Euler => IntegratorOpts::Euler(
            parse_euler_opts(&raw_int_opts).unwrap(),
        ),
        IntegratorType::EulerDebug => IntegratorOpts::EulerDebug(
            parse_euler_opts(&raw_int_opts).unwrap(),
        ),
        IntegratorType::Rkdp5 => IntegratorOpts::Rkdp5(
            parse_rkdp5_opts(&raw_int_opts).unwrap(),
        ),
    }
}

#[derive(Clone, Debug)]
pub struct ExperimentArgs {
    pub toml_name: String,
    pub ty: ExperimentType,
    pub final_t: Time,
    pub cil_mag: f64,
    pub coa_mag: Option<f64>,
    pub cal_mag: Option<f64>,
    pub adh_scale: Option<f64>,
    pub snap_period: Time,
    pub max_on_ram: usize,
    pub seeds: Vec<Option<u64>>,
    pub int_opts: IntegratorOpts,
}

impl TryFrom<&PathBuf> for ExperimentArgs {
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
        let raw_exp_args = {
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
        let ty = get_exp_type(&raw_exp_args)?;
        let final_t = Time({
            let value = get_value(&raw_exp_args, "final_t")?;
            parse_f64(&value)?
        });
        let cil = {
            let value = get_value(&raw_exp_args, "cil")?;
            parse_f64(&value)?
        };
        let coa = get_optional_f64(&raw_exp_args, "coa")?;
        let adh = get_optional_f64(&raw_exp_args, "adh")?;
        let cal = get_optional_f64(&raw_exp_args, "cal")?;
        let snap_period = Time(
            get_optional_f64(&raw_exp_args, "snap_period")?
                .unwrap_or(5.0),
        );
        let max_on_ram =
            get_optional_usize(&raw_exp_args, "max_on_ram")?
                .unwrap_or(1000);
        let seeds = get_seeds(&raw_exp_args)?;
        let toml_name: String = toml_path
            .file_stem()
            .unwrap_or_else(|| {
                panic!(
                    "could not extract file name from path: {}",
                    &toml_path.to_str().unwrap()
                )
            })
            .to_str()
            .unwrap_or_else(|| {
                panic!(
                    "could not convert file name to string from path: {}",
                    &toml_path.to_str().unwrap()
                )
            })
            .into();
        let int_opts: IntegratorOpts =
            get_integrator_opts(&raw_exp_args);
        Ok(ExperimentArgs {
            toml_name,
            final_t,
            cil_mag: cil,
            coa_mag: coa,
            cal_mag: cal,
            adh_scale: adh,
            snap_period,
            max_on_ram,
            seeds,
            int_opts,
            ty,
        })
    }
}
