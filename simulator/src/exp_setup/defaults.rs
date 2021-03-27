use crate::cell::chemistry::{distrib_gens, RgtpDistribution};
use crate::exp_setup::markers::ALL;
use crate::parameters::quantity::{
    Force, General, Length, Quantity, Stress, Time, Tinv, Viscosity,
};
use crate::parameters::{
    CharQuantities, RawCoaParams, RawInteractionParams,
    RawParameters, RawPhysicalContactParams, RawWorldParameters,
};
use crate::NVERTS;
use once_cell::sync::Lazy;

//TODO: Document all justifications for characteristic quantities, rather than
// having reader refer to the SI.
pub static CHAR_FORCE: Lazy<Force> = Lazy::new(|| {
    (Stress(1.0).kilo().g()
        * Length(100.0).nano().g()
        * Length(10.0).micro().g())
    .to_force()
    .unwrap()
});
pub const CHAR_VISCOSITY: Viscosity = Viscosity(0.1);
pub static CHAR_LENGTH: Lazy<Length> =
    Lazy::new(|| Length(1.0).micro());
pub const CHAR_TIME: Time = Time(2.0);
pub static CHAR_L3D: Lazy<Length> =
    Lazy::new(|| Length(10.0).micro());
pub const CHAR_KGTP: Tinv = Tinv(1e-4);

/// Default characteristic quantities.Refer to SI of first two papers for
/// justification of the values used.
pub static CHAR_QUANTS: Lazy<CharQuantities> =
    Lazy::new(|| CharQuantities {
        eta: CHAR_VISCOSITY,
        l: *CHAR_LENGTH,
        t: CHAR_TIME,
        f: *CHAR_FORCE,
        l3d: *CHAR_L3D,
        kgtp: CHAR_KGTP,
    });

pub static MAX_CELL_V: Lazy<General> =
    Lazy::new(|| Length(3.0).micro().g() * Tinv(1.0 / 60.0).g());
pub const ADH_INDEX: f64 = 0.99;
pub static ADH_MAG: Lazy<Force> = Lazy::new(|| {
    (CHAR_VISCOSITY.g() * (*MAX_CELL_V))
        .scale(1.0 / NVERTS as f64)
        .to_force()
        .expect(
            "Procedure for generating default force does \
             not produce a force. Check units!",
        )
});
pub static CELL_DIAMETER: Lazy<Length> =
    Lazy::new(|| Length(40.0).micro());

/// Default raw parameters for cells.
pub static RAW_PARAMS: Lazy<RawParameters> = Lazy::new(|| {
    let rgtp_d = (Length(0.1_f64.sqrt()).micro().pow(2.0).g()
        / Time(1.0).g())
    .to_diffusion()
    .unwrap();
    let init_rac = RgtpDistribution::new(
        distrib_gens::specific_uniform(0.1, ALL),
        distrib_gens::specific_uniform(0.1, ALL),
    );
    RawParameters {
        cell_diam: *CELL_DIAMETER,
        stiffness_cortex: Stress(8.0).kilo(),
        lm_h: Length(200.0).nano(),
        halfmax_rgtp_max_f_frac: 0.3,
        halfmax_rgtp_frac: 0.4,
        lm_ss: Stress(10.0).kilo(),
        rho_friction: 0.2,
        stiffness_cyto: Force(1e-7),
        diffusion_rgtp: rgtp_d,
        k_mem_off: Tinv(0.15),
        k_mem_on: Tinv(0.02),
        kgtp_rac: Tinv(1e-4).scale(24.0),
        kgtp_rac_auto: Tinv(1e-4).scale(500.0),
        kdgtp_rac: Tinv(1e-4).scale(8.0),
        kdgtp_rho_on_rac: Tinv(1e-4).scale(4000.0),
        halfmax_tension_inhib: 0.1,
        tension_inhib: 40.0,
        kgtp_rho: Tinv(1e-4).scale(28.0),
        kgtp_auto_rho: Tinv(1e-4).scale(390.0),
        kdgtp_rho: Tinv(1e-4).scale(60.0),
        kdgtp_rac_on_rho: Tinv(1e-4).scale(400.0),
        randomization: false,
        rand_avg_t: Time(40.0 * 60.0),
        rand_std_t: Time(0.1 * 40.0 * 60.0),
        rand_mag: 10.0,
        rand_vs: 0.25,
        init_rac,
        init_rho: init_rac,
    }
});

/// We take the viscosity of the world to be 0.29 N m^-2. We
/// divide viscosity by the number of vertices, on a cell in
/// order to scale it properly.
pub fn vertex_viscosity(char_quants: &CharQuantities) -> Viscosity {
    char_quants.eta.scale(2.9 / (NVERTS as f64))
}

pub static PHYS_CLOSE_DIST: Lazy<Length> =
    Lazy::new(|| Length(0.5).micro());
pub static PHYS_CLOSE_DIST_ONE_AT: Lazy<Length> =
    Lazy::new(|| *PHYS_CLOSE_DIST);
pub static PHYS_CLOSE_DIST_ZERO_AT: Lazy<Length> =
    Lazy::new(|| PHYS_CLOSE_DIST.scale(3.0));

pub const CIL_MAG: f64 = 60.0;

pub const COA_LOS_PENALTY: f64 = 2.0;
pub static COA_HALFMAX_DIST: Lazy<Length> =
    Lazy::new(|| Length(110.0).micro());

pub static RAW_COA_PARAMS_WITH_ZERO_MAG: Lazy<RawCoaParams> =
    Lazy::new(|| RawCoaParams {
        los_penalty: COA_LOS_PENALTY,
        halfmax_dist: *COA_HALFMAX_DIST,
        mag: 0.0,
        too_close_dist: Length(1.0).micro(), //PHYS_CLOSE_DIST.scale(2.0),
    });

pub static RAW_WORLD_PARAMS: Lazy<RawWorldParameters> =
    Lazy::new(|| {
        let one_at = *PHYS_CLOSE_DIST;
        RawWorldParameters {
            vertex_eta: vertex_viscosity(&CHAR_QUANTS),
            interactions: RawInteractionParams {
                coa: None,
                chem_attr: None,
                bdry: None,
                phys_contact: RawPhysicalContactParams {
                    crl_one_at: one_at,
                    zero_at: one_at.scale(2.0),
                    adh_mag: None,
                    cal_mag: None,
                    cil_mag: CIL_MAG,
                    adh_break: None,
                },
            },
        }
    });
