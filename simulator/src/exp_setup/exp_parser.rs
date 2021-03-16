use crate::cell::chemistry::distrib_gens::{
    random, specific_random, specific_uniform,
};
use crate::exp_setup::markers::mark_verts;
use crate::exp_setup::ExperimentType;
use crate::parameters::quantity::Time;
use crate::toml_parse::{
    get_optional_f64, get_optional_path, get_optional_usize,
    get_value, parse_array, parse_bool, parse_f64, parse_string,
    parse_table, parse_u64, ParseErr,
};
use crate::utils::pcg32::Pcg32;
use crate::world::{EulerOpts, IntegratorOpts, Rkdp5Opts};
use crate::NVERTS;
use serde::{Deserialize, Serialize};
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
            rgtp_distrib_defs_per_cell: {
                let v = get_value(
                    raw_exp_args,
                    "rgtp_distrib_defs_per_cell",
                )?;
                let arr = parse_array(&v)?;
                let mut r = vec![];
                for value in arr.iter() {
                    r.push(parse_rgtp_distrib_def(&value)?);
                }
                r
            },
        },
        "PyCompare" => ExperimentType::PyCompare {
            num_cells: {
                let v = get_value(raw_exp_args, "num_cells")?;
                (parse_u64(&v)?) as usize
            },
            py_main: get_optional_path(raw_exp_args, "py_main")?,
        },
        "NCells" => ExperimentType::NCells {
            num_cells: {
                let v = get_value(raw_exp_args, "num_cells")?;
                (parse_u64(&v)?) as usize
            },
        },
        _ => {
            return Err(ParseErr::UnknownExperiment(exp_type_str));
        }
    };
    Ok(exp_type)
}

fn get_randomization(raw_exp_args: &Table) -> Result<bool, ParseErr> {
    let value = get_value(raw_exp_args, "randomization")?;
    parse_bool(&value)
}

fn get_seeds(raw_exp_args: &Table) -> Result<Vec<u64>, ParseErr> {
    let value = get_value(raw_exp_args, "seeds")?;
    let raw_seeds = parse_array(&value)?;
    let mut seeds: Vec<u64> = Vec::with_capacity(raw_seeds.len());
    for raw_seed in raw_seeds {
        let seed = parse_u64(&raw_seed)?;
        seeds.push(seed);
    }
    if seeds.len() > 0 {
        Ok(seeds)
    } else {
        Err(ParseErr::NoSeed)
    }
}

fn get_marked_verts(
    raw_distrib_defn: &Table,
) -> Result<Vec<usize>, ParseErr> {
    let value = get_value(raw_distrib_defn, "marked_verts")?;
    let raw_verts = parse_array(&value)?;
    if raw_verts.len() > NVERTS {
        Err(ParseErr::TooManyMarkedVerts(NVERTS, raw_verts.len()))
    } else {
        let mut verts = Vec::with_capacity(raw_verts.len());
        for raw_vert in raw_verts {
            let vert = parse_u64(&raw_vert)? as usize;
            verts.push(vert)
        }
        Ok(verts)
    }
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

pub fn get_rgtp_distrib(
    raw_exp_args: &Table,
    distrib_name: &str,
) -> DistribDef {
    let distrib_table_key = format!("init_{}", distrib_name);
    get_value(raw_exp_args, &distrib_table_key).map_or_else(
        |_| DistribDef::default(),
        |value| {
            let raw_distrib = parse_table(&value).unwrap();
            let ty_value = get_value(&raw_distrib, "ty").unwrap();
            let ty = parse_string(&ty_value).unwrap();
            let frac = get_value(&raw_distrib, "frac")
                .map(|inner| parse_f64(&inner).unwrap())
                .unwrap();
            match ty.as_str() {
                "Random" => DistribDef::Random { frac },
                "SpecificRandom" => DistribDef::SpecificRandom {
                    frac,
                    marked_verts: get_marked_verts(&raw_distrib)
                        .unwrap(),
                },
                "SpecificUniform" => {
                    DistribDef::SpecificUniform {
                        frac,
                        marked_verts: get_marked_verts(&raw_distrib)
                            .unwrap(),
                    }
                }
                s => Err(ParseErr::UnknownRgtpDistribDefn(s.into()))
                    .unwrap(),
            }
        },
    )
}

pub fn get_rgtp_distribs(raw_exp_args: &Table) -> RgtpDistribDefs {
    let rac_acts = get_rgtp_distrib(&raw_exp_args, "rac_acts");
    let rac_inacts = get_rgtp_distrib(&raw_exp_args, "rac_inacts");
    let rho_acts = get_rgtp_distrib(&raw_exp_args, "rho_acts");
    let rho_inacts = get_rgtp_distrib(&raw_exp_args, "rho_inacts");

    RgtpDistribDefs {
        rac_acts,
        rac_inacts,
        rho_acts,
        rho_inacts,
    }
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
    pub randomization: bool,
    pub seeds: Vec<u64>,
    pub int_opts: IntegratorOpts,
    pub rgtp_distrib_defs: RgtpDistribDefs,
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
        let randomization = get_randomization(&raw_exp_args)?;
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
        let rgtp_distribs = get_rgtp_distribs(&raw_exp_args);
        Ok(ExperimentArgs {
            toml_name,
            final_t,
            cil_mag: cil,
            coa_mag: coa,
            cal_mag: cal,
            adh_scale: adh,
            snap_period,
            max_on_ram,
            rgtp_distrib_defs: rgtp_distribs,
            randomization,
            seeds,
            int_opts,
            ty,
        })
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum DistribDef {
    Random { frac: f64 },
    SpecificRandom { frac: f64, marked_verts: Vec<usize> },
    SpecificUniform { frac: f64, marked_verts: Vec<usize> },
}

impl DistribDef {
    pub fn into_distrib(self, rng: &mut Pcg32) -> [f64; NVERTS] {
        match self {
            DistribDef::Random { frac } => random(rng, frac),
            DistribDef::SpecificRandom {
                frac,
                marked_verts,
            } => specific_random(rng, frac, mark_verts(marked_verts)),
            DistribDef::SpecificUniform {
                frac,
                marked_verts,
            } => specific_uniform(frac, mark_verts(marked_verts)),
        }
    }
}

impl Default for DistribDef {
    fn default() -> Self {
        DistribDef::Random { frac: 0.1 }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct RgtpDistribDef {
    acts: DistribDef,
    inacts: DistribDef,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct RgtpDistribDefs {
    rac: RgtpDistribDef,
    rho: RgtpDistribDef,
}

