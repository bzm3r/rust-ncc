use crate::cell::chemistry::{
    calc_conc_rgtps, calc_kdgtps_rac, calc_kdgtps_rho,
    calc_kgtps_rac, calc_kgtps_rho, calc_net_fluxes, RacRandState,
    RgtpDistribution,
};
use crate::cell::mechanics::{
    calc_cyto_forces, calc_edge_forces, calc_edge_vecs,
    calc_rgtp_forces,
};
use crate::interactions::{
    Contact, Interactions, RelativeRgtpActivity,
};
use crate::math::geometry::{
    is_point_in_poly, lsegs_intersect, LineSeg2D, Poly,
};
use crate::math::v2d::{SqP2d, V2d};
use crate::math::{hill_function3, max_f64};
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

/// `CoreState` contains all the variables that are simulated between geometric
/// updates. They are simulated using ODEs which are then integrated using
/// either the Euler method or Runge-Kutta Dormand-Prince 5 (Matlab's `ode45`).
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct Core {
    /// Polygon representing cell shape.
    pub poly: [V2d; NVERTS],
    /// Fraction of Rac1 active at each vertex.
    pub rac_acts: [f64; NVERTS],
    /// Fraction of Rac1 inactive at each vertex.
    pub rac_inacts: [f64; NVERTS],
    /// Fraction of RhoA active at each vertex.
    pub rho_acts: [f64; NVERTS],
    /// Fraction of RhoA inactive at each vertex.
    pub rho_inacts: [f64; NVERTS],
    /// Geometric state resulting from this core state.
    pub geom: GeomState,
}

/// `DCoreDt` is the derivative of `CoreState`.
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct DCoreDt {
    /// Polygon representing cell shape.
    pub poly: [V2d; NVERTS],
    /// Fraction of Rac1 active at each vertex.
    pub rac_acts: [f64; NVERTS],
    /// Fraction of Rac1 inactive at each vertex.
    pub rac_inacts: [f64; NVERTS],
    /// Fraction of RhoA active at each vertex.
    pub rho_acts: [f64; NVERTS],
    /// Fraction of RhoA inactive at each vertex.
    pub rho_inacts: [f64; NVERTS],
}

impl From<&Core> for DCoreDt {
    fn from(_: &Core) -> Self {
        unimplemented!()
    }
}

impl DCoreDt {
    pub fn time_step(&self, dt: f64) -> Core {
        let mut poly = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = dt * self.poly[i];
            rac_acts[i] = dt * self.rac_acts[i];
            rac_inacts[i] = dt * self.rac_inacts[i];
            rho_acts[i] = dt * self.rho_acts[i];
            rho_inacts[i] = dt * self.rho_inacts[i];
        }

        Core::new(poly, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }
}

/// `PowCoreState` results from multiplication/division of `CoreState`s, or
/// taking their power using `CoreState::powi`.
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct SqCore {
    /// Polygon representing cell shape.
    pub poly: [SqP2d; NVERTS],
    /// Fraction of Rac1 active at each vertex.
    pub rac_acts: [f64; NVERTS],
    /// Fraction of Rac1 inactive at each vertex.
    pub rac_inacts: [f64; NVERTS],
    /// Fraction of RhoA active at each vertex.
    pub rho_acts: [f64; NVERTS],
    /// Fraction of RhoA inactive at each vertex.
    pub rho_inacts: [f64; NVERTS],
}

impl Add for Core {
    type Output = Core;

    fn add(self, rhs: Core) -> Core {
        let mut poly = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self.poly[i] + rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] + rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] + rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] + rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] + rhs.rho_inacts[i]
        }

        Core::new(poly, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }
}

impl Mul for Core {
    type Output = SqCore;

    fn mul(self, rhs: Core) -> Self::Output {
        let mut poly = [SqP2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self.poly[i] * rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] * rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] * rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] * rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] * rhs.rho_inacts[i];
        }

        Self::Output {
            poly,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Mul<&Core> for &Core {
    type Output = SqCore;

    fn mul(self, rhs: &Core) -> Self::Output {
        let mut poly = [SqP2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self.poly[i] * rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] * rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] * rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] * rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] * rhs.rho_inacts[i];
        }

        Self::Output {
            poly,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Sub for Core {
    type Output = Core;

    fn sub(self, rhs: Core) -> Core {
        let mut poly = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self.poly[i] - rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] - rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] - rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] - rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] - rhs.rho_inacts[i]
        }

        Core::new(poly, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }
}

impl Mul<f64> for Core {
    type Output = Core;

    fn mul(self, rhs: f64) -> Self::Output {
        let mut poly = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self.poly[i].scale(rhs);
            rac_acts[i] = self.rac_acts[i] * rhs;
            rac_inacts[i] = self.rac_inacts[i] * rhs;
            rho_acts[i] = self.rho_acts[i] * rhs;
            rho_inacts[i] = self.rho_inacts[i] * rhs;
        }

        Core::new(poly, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }
}

impl Mul<Core> for f64 {
    type Output = Core;

    fn mul(self, rhs: Core) -> Self::Output {
        rhs * self
    }
}

impl Add<f64> for Core {
    type Output = Core;

    fn add(self, rhs: f64) -> Self::Output {
        let mut poly = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self.poly[i] + rhs;
            rac_acts[i] = self.rac_acts[i] + rhs;
            rac_inacts[i] = self.rac_inacts[i] + rhs;
            rho_acts[i] = self.rho_acts[i] + rhs;
            rho_inacts[i] = self.rho_inacts[i] + rhs;
        }

        Core::new(poly, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }
}

impl Add<Core> for f64 {
    type Output = Core;

    fn add(self, rhs: Core) -> Self::Output {
        rhs + self
    }
}

impl Div<Core> for SqCore {
    type Output = Core;

    fn div(self, rhs: Core) -> Self::Output {
        let mut poly = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self.poly[i] / rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] / rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] / rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] / rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] / rhs.rho_inacts[i]
        }

        Core::new(poly, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }
}

/// Records the mechanical state of a cell.
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct MechState {
    /// Strain each edge is under, where resting edge length is
    /// defined in the cell's parameters.
    pub edge_strains: [f64; NVERTS],
    /// Forces on each vertex due to Rho GTPase activity.
    pub rgtp_forces: [V2d; NVERTS],
    /// Forces on each vertex due to cytoplasmic pressure.
    pub cyto_forces: [V2d; NVERTS],
    /// Forces on each vertex due to edge-edge (elastic) forces.
    pub edge_forces: [V2d; NVERTS],
    /// Average of the strain in edges which are under tension (i.e.
    /// they are longer than their initial resting edge length.
    pub avg_tens_strain: f64,
    /// Sum of all forces that are acting on a vertex, except for
    /// adhesion, which comes from interaction information.
    pub sum_forces: [V2d; NVERTS],
}

/// Calculates the various rates necessary to define the ODEs
/// simulating biochemistry. These are (`X` is either `rac` or `rho`):
///     * `kdgtp_X`: Rho GTPase inactivation rates
///     * `kgtp_X`: Rho GTPase activation rates
///     * `X_act_net_fluxes`: diffusion fluxes between vertices of
/// active form of Rho GTPase
///     * `X_inact_net_fluxes`: diffusion fluxes between vertices of
/// inactive form of Rho GTPase
///     * `X_cyto`: fraction of Rho GTPase in the cytoplasm
///     * `x_tens`: "tension" factor that affects Rac1 activation
/// rate, calculated based on average tensile strain in cell (i.e.
/// how stretched the cell is).
#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct ChemState {
    pub kdgtps_rac: [f64; NVERTS],
    pub kgtps_rac: [f64; NVERTS],
    pub rac_act_net_fluxes: [f64; NVERTS],
    pub rac_inact_net_fluxes: [f64; NVERTS],
    pub kdgtps_rho: [f64; NVERTS],
    pub kgtps_rho: [f64; NVERTS],
    pub rac_cyto: f64,
    pub rho_cyto: f64,
    pub rho_act_net_fluxes: [f64; NVERTS],
    pub rho_inact_net_fluxes: [f64; NVERTS],
    pub x_tens: f64,
}

#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct GeomState {
    /// Unit edge vectors which point from position of vertex `vi`
    /// to position of vertex `vi + 1`, where `vi + 1` is calculated
    /// modulo `NVERTS`. Note that an edge is defined by its "lower"
    /// index, modulo `NVERTS`. That is, the edge `(0, 1)`, is
    /// different from the edge `(1, 2)`, and `(0, 1)` is also
    /// different from `(15, 0)` (assuming that `NVERTS == 16` in
    /// this example).
    pub unit_edge_vecs: [V2d; NVERTS],
    /// Length of edges. Each edge is defined by its smallest vertex.
    pub edge_lens: [f64; NVERTS],
    /// Inward pointing unit vectors at each vertex. These are
    /// calculated so that they bisect the angle between the two
    /// edges which meet at a vertex.
    pub unit_in_vecs: [V2d; NVERTS],
}

impl From<&[V2d; NVERTS]> for GeomState {
    fn from(poly: &[V2d; 16]) -> Self {
        // Calculate edge vectors of a polygon.
        let evs = calc_edge_vecs(poly);
        // Calculate magnitude of each edge vec, to get its length.
        let mut edge_lens = [0.0_f64; NVERTS];
        (0..NVERTS).for_each(|i| edge_lens[i] = (&evs[i]).mag());
        // Divide each edge vector by its magnitude to get the
        // corresponding unit vector.
        let mut unit_edge_vecs = [V2d::default(); NVERTS];
        (0..NVERTS)
            .for_each(|i| unit_edge_vecs[i] = (&evs[i]).unitize());
        // Given two unit edge vectors, find the vector which points
        // into the polygon and bisects the angle
        let mut unit_in_vecs = [V2d::default(); NVERTS];
        (0..NVERTS).for_each(|i| {
            let im1 = circ_ix_minus(i, NVERTS);
            let tangent =
                (unit_edge_vecs[i] + unit_edge_vecs[im1]).unitize();
            unit_in_vecs[i] = tangent.normal();
        });

        GeomState {
            unit_edge_vecs,
            edge_lens,
            unit_in_vecs,
        }
    }
}

pub fn fmt_var_arr<T: fmt::Display>(
    f: &mut fmt::Formatter<'_>,
    description: &str,
    vars: &[T; NVERTS],
) -> fmt::Result {
    let contents = vars
        .iter()
        .map(|x| format!("{}", x))
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(f, "{}: [{}]", description, contents)
}

impl Display for Core {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_var_arr(f, "vertex_coords", &self.poly)?;
        fmt_var_arr(f, "rac_acts", &self.rac_acts)?;
        fmt_var_arr(f, "rac_inacts", &self.rac_inacts)?;
        fmt_var_arr(f, "rho_acts", &self.rho_acts)?;
        fmt_var_arr(f, "rho_inacts", &self.rho_inacts)
    }
}

impl Display for GeomState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_var_arr(f, "edge_lens", &self.edge_lens)?;
        fmt_var_arr(f, "uivs", &self.unit_in_vecs)?;
        fmt_var_arr(f, "uevs", &self.unit_edge_vecs)
    }
}

impl Display for MechState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_var_arr(f, "rgtp_forces", &self.rgtp_forces)?;
        fmt_var_arr(f, "edge_strains", &self.edge_strains)?;
        writeln!(f, "avg_tens_strain: {}", self.avg_tens_strain)?;
        fmt_var_arr(f, "edge_forces", &self.edge_forces)?;
        fmt_var_arr(f, "cyto_forces", &self.cyto_forces)?;
        fmt_var_arr(f, "tot_forces", &self.sum_forces)
    }
}

impl Display for ChemState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "x_tens: {}", self.x_tens)?;
        fmt_var_arr(f, "kgtps_rac", &self.kgtps_rac)?;
        fmt_var_arr(f, "kdgtps_rac", &self.kdgtps_rac)?;
        fmt_var_arr(f, "kgtps_rho", &self.kgtps_rho)?;
        fmt_var_arr(f, "kdgtps_rho", &self.kdgtps_rac)
    }
}

impl Core {
    pub fn calc_mech_state(
        &self,
        parameters: &Parameters,
    ) -> MechState {
        let GeomState {
            unit_edge_vecs,
            edge_lens,
            unit_in_vecs,
        } = &self.geom;
        let rgtp_forces = calc_rgtp_forces(
            &self.rac_acts,
            &self.rho_acts,
            unit_in_vecs,
            parameters.halfmax_vertex_rgtp,
            parameters.const_protrusive,
            parameters.const_retractive,
        );
        let cyto_forces = calc_cyto_forces(
            &self.poly,
            unit_in_vecs,
            parameters.rest_area,
            parameters.stiffness_cyto,
        );
        // Calculate strain in each edge.
        let mut edge_strains = [0.0_f64; NVERTS];
        (0..NVERTS).for_each(|i| {
            edge_strains[i] =
                (edge_lens[i] / parameters.rest_edge_len) - 1.0
        });
        let edge_forces = calc_edge_forces(
            &edge_strains,
            unit_edge_vecs,
            parameters.stiffness_edge,
        );
        // If strain is positive (tensile), then consider it in this
        // averaging, otherwise don't. This is because I'm assuming
        // that compression does not have an effect on Rac1 activity.
        // Only tension is considered to have an effect (see refs. in
        // SI.
        //TODO(BM): what is the latest on this front? Recent paper
        // (ELife?) which suggests not true for migrating cells.
        let avg_tens_strain = edge_strains
            .iter()
            .map(|&es| if es < 0.0 { 0.0 } else { es })
            .sum::<f64>()
            / NVERTS as f64;
        // Sum of all the non-adhesive forces acting on the cell.
        let mut sum_forces = [V2d::default(); NVERTS];
        (0..NVERTS).for_each(|i| {
            sum_forces[i] =
                rgtp_forces[i] + cyto_forces[i] + edge_forces[i]
                    - edge_forces[circ_ix_minus(i, NVERTS)];
        });
        MechState {
            edge_strains,
            rgtp_forces,
            cyto_forces,
            edge_forces,
            avg_tens_strain,
            sum_forces,
        }
    }

    pub fn calc_chem_state(
        &self,
        mech_state: &MechState,
        rac_rand_state: &RacRandState,
        interactions: &Interactions,
        parameters: &Parameters,
    ) -> ChemState {
        let GeomState { edge_lens, .. } = self.geom;
        // Need to calculate average length of edges meeting at
        // a vertex in order to roughly approximate diffusion related
        // flux of Rho GTPase from neighbouring vertices. Provides
        // an approximation for the length of the membrane abstracted
        // by the edges that meet at that vertex.
        let mut avg_edge_lens: [f64; NVERTS] = [0.0_f64; NVERTS];
        (0..NVERTS).for_each(|i| {
            let im1 = circ_ix_minus(i, NVERTS);
            avg_edge_lens[i] = (edge_lens[i] + edge_lens[im1]) / 2.0;
        });

        // Concentration of Rho GTPase (active/inactive), used to
        // calculate diffusive flux between vertices.
        let conc_rac_acts =
            calc_conc_rgtps(&avg_edge_lens, &self.rac_acts);
        let conc_rac_inacts =
            calc_conc_rgtps(&avg_edge_lens, &self.rac_inacts);
        let conc_rho_acts =
            calc_conc_rgtps(&avg_edge_lens, &self.rho_acts);
        let conc_rho_inacts =
            calc_conc_rgtps(&avg_edge_lens, &self.rho_inacts);

        let kgtps_rac = calc_kgtps_rac(
            &self.rac_acts,
            &conc_rac_acts,
            &rac_rand_state.x_rands,
            &interactions.x_coas,
            &interactions.x_cils,
            &interactions.x_chem_attrs,
            &interactions.x_cals,
            parameters.kgtp_rac,
            parameters.kgtp_rac_auto,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let MechState {
            avg_tens_strain, ..
        } = mech_state;
        let x_tens = parameters.tension_inhib
            * hill_function3(
                parameters.halfmax_tension_inhib,
                *avg_tens_strain,
            );
        let kdgtps_rac = calc_kdgtps_rac(
            &self.rac_acts,
            &conc_rho_acts,
            &interactions.x_cils,
            x_tens,
            parameters.kdgtp_rac,
            parameters.kdgtp_rho_on_rac,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let kgtps_rho = calc_kgtps_rho(
            &self.rho_acts,
            &conc_rho_acts,
            &interactions.x_cils,
            parameters.kgtp_rho,
            parameters.halfmax_vertex_rgtp_conc,
            parameters.kgtp_rho_auto,
        );
        let kdgtps_rho = calc_kdgtps_rho(
            &self.rho_acts,
            &conc_rac_acts,
            parameters.kdgtp_rho,
            parameters.kdgtp_rac_on_rho,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let rac_act_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rac_acts,
        );
        let rho_act_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rho_acts,
        );
        let rac_inact_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rac_inacts,
        );
        let rho_inact_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rho_inacts,
        );

        let rac_cyto = 1.0
            - self.rac_acts.iter().sum::<f64>()
            - self.rac_inacts.iter().sum::<f64>();
        let rho_cyto = 1.0
            - self.rho_acts.iter().sum::<f64>()
            - self.rho_inacts.iter().sum::<f64>();
        ChemState {
            x_tens,
            kdgtps_rac,
            kgtps_rac,
            rac_act_net_fluxes,
            rac_inact_net_fluxes,
            rho_act_net_fluxes,
            rho_inact_net_fluxes,
            kdgtps_rho,
            kgtps_rho,
            rac_cyto,
            rho_cyto,
        }
    }

    /// Calculate the right hand side of the ODEs simulating cell
    /// vertex motion and biochemistry.
    pub fn derivative(
        &self,
        rac_rand_state: &RacRandState,
        interactions: &Interactions,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> DCoreDt {
        //TODO: is it necessary to recalculate chem/mech/geom states in
        // `derivative`, if we have saved this info in the `Cell` struct?
        // What is the importance of `interactions`---might it have changed
        // since the last time we calculated these?
        let mech_state = self.calc_mech_state(parameters);
        let chem_state = self.calc_chem_state(
            &mech_state,
            rac_rand_state,
            interactions,
            parameters,
        );
        let mut delta = DCoreDt::default();
        for i in 0..NVERTS {
            // rate of rac deactivation * current fraction of rac active
            let inactivated_rac =
                chem_state.kdgtps_rac[i] * self.rac_acts[i];
            let activated_rac =
                chem_state.kgtps_rac[i] * self.rac_inacts[i];
            let inactivated_rho =
                chem_state.kdgtps_rho[i] * self.rho_acts[i];
            let activated_rho =
                chem_state.kgtps_rho[i] * self.rho_inacts[i];

            let delta_rac_activated = activated_rac - inactivated_rac;
            let delta_rho_activated = activated_rho - inactivated_rho;

            let rac_cyto_exchange = {
                let rac_mem_on =
                    parameters.k_mem_on_vertex * chem_state.rac_cyto;
                let rac_mem_off =
                    parameters.k_mem_off * self.rac_inacts[i];
                rac_mem_on - rac_mem_off
            };
            let rho_cyto_exchange = {
                let rho_mem_on =
                    parameters.k_mem_on_vertex * chem_state.rho_cyto;
                let rho_mem_off =
                    parameters.k_mem_off * self.rho_inacts[i];
                rho_mem_on - rho_mem_off
            };

            let vertex_rac_act_flux =
                chem_state.rac_act_net_fluxes[i];
            let vertex_rac_inact_flux =
                chem_state.rac_inact_net_fluxes[i];
            let vertex_rho_act_flux =
                chem_state.rho_act_net_fluxes[i];
            let vertex_rho_inact_flux =
                chem_state.rho_inact_net_fluxes[i];

            delta.rac_acts[i] =
                delta_rac_activated + vertex_rac_act_flux;
            delta.rac_inacts[i] = rac_cyto_exchange
                + vertex_rac_inact_flux
                - delta_rac_activated;
            delta.rho_acts[i] =
                delta_rho_activated + vertex_rho_act_flux;
            delta.rho_inacts[i] = rho_cyto_exchange
                + vertex_rho_inact_flux
                - delta_rho_activated;
            delta.poly[i] = (1.0 / world_parameters.vertex_eta)
                * (mech_state.sum_forces[i] + interactions.x_adhs[i]);
        }
        delta
    }

    pub fn init(
        poly: [V2d; NVERTS],
        init_rac: RgtpDistribution,
        init_rho: RgtpDistribution,
    ) -> Core {
        Core::new(
            poly,
            init_rac.active,
            init_rac.inactive,
            init_rho.active,
            init_rho.inactive,
        )
    }

    pub fn new(
        poly: [V2d; NVERTS],
        rac_acts: [f64; NVERTS],
        rac_inacts: [f64; NVERTS],
        rho_acts: [f64; NVERTS],
        rho_inacts: [f64; NVERTS],
    ) -> Core {
        let geom = GeomState::from(&poly);
        Core {
            poly,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
            geom,
        }
    }

    pub fn abs(&self) -> Core {
        let mut vertex_coords = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].abs();
            rac_acts[i] = self.rac_acts[i].abs();
            rac_inacts[i] = self.rac_inacts[i].abs();
            rho_acts[i] = self.rho_acts[i].abs();
            rho_inacts[i] = self.rho_inacts[i].abs();
        }

        Core::new(
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        )
    }

    pub fn square(&self) -> SqCore {
        self * self
    }

    pub fn max(&self, other: &Core) -> Core {
        let mut poly = [V2d::default(); NVERTS];
        let mut rac_acts = [0.0_f64; NVERTS];
        let mut rac_inacts = [0.0_f64; NVERTS];
        let mut rho_acts = [0.0_f64; NVERTS];
        let mut rho_inacts = [0.0_f64; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = poly[i].max(&other.poly[i]);
            rac_acts[i] =
                max_f64(self.rac_acts[i], other.rac_acts[i]);
            rac_inacts[i] =
                max_f64(self.rac_inacts[i], other.rac_inacts[i]);
            rho_acts[i] =
                max_f64(self.rho_acts[i], other.rho_acts[i]);
            rho_inacts[i] =
                max_f64(self.rho_inacts[i], other.rho_inacts[i]);
        }

        Core::new(poly, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }

    /// Calculate which Rho GTPase has dominates in terms of effect
    /// at this vertex.
    pub fn calc_relative_rgtp_activity(
        &self,
        parameters: &Parameters,
    ) -> [RelativeRgtpActivity; NVERTS] {
        let mut r = [RelativeRgtpActivity::RhoDominant(0.0); NVERTS];
        self.rac_acts
            .iter()
            .zip(self.rho_acts.iter())
            .enumerate()
            .for_each(|(ix, (&rac, &rho))| {
                r[ix] = RelativeRgtpActivity::from_f64(
                    hill_function3(
                        parameters.halfmax_vertex_rgtp,
                        rac,
                    ) - hill_function3(
                        parameters.halfmax_vertex_rgtp,
                        rho,
                    ),
                );
            });
        r
    }

    #[cfg(feature = "validate")]
    pub fn validate(&self, loc_str: &str) -> Result<(), String> {
        if self.rac_acts.iter().any(|&r| r < 0.0_f64) {
            return Err(format!(
                "{}: neg rac_acts: {:?}",
                loc_str, self.rac_acts
            ));
        }
        if self.rac_inacts.iter().any(|&r| r < 0.0_f64) {
            return Err(format!(
                "{}: neg rac_inacts: {:?}",
                loc_str, self.rac_inacts
            ));
        }
        if self.rho_acts.iter().any(|&r| r < 0.0_f64) {
            return Err(format!(
                "{}: neg rho_acts: {:?}",
                loc_str, self.rho_acts
            ));
        }
        if self.rho_inacts.iter().any(|&r| r < 0.0_f64) {
            return Err(format!(
                "{}: neg rho_inacts: {:?}",
                loc_str, self.rho_inacts
            ));
        }
        let sum_rac_mem = self.rac_inacts.iter().sum::<f64>()
            + self.rac_acts.iter().sum::<f64>();
        if !(0.0..=1.0).contains(&sum_rac_mem) {
            return Err(format!(
                "{}: problem in sum of rac_mem: {}",
                loc_str, sum_rac_mem
            ));
        }
        let sum_rho_mem = self.rho_inacts.iter().sum::<f64>()
            + self.rho_acts.iter().sum::<f64>();
        if !(0.0..=1.0).contains(&sum_rho_mem) {
            return Err(format!(
                "{}: problem in sum of rho_mem: {}",
                loc_str, sum_rho_mem
            ));
        }
        Ok(())
        // println!("{}: successfully validated", loc_str)
    }

    pub fn strict_enforce_volume_exclusion(
        &mut self,
        old_vs: &[V2d; NVERTS],
        contacts: &[Contact],
    ) -> Result<(), VolExErr> {
        // confirm_volume_exclusion(&old_vs, &contacts, "old_vs")?;

        self.enforce_volume_exclusion(old_vs, contacts)?;

        // confirm_volume_exclusion(&self.poly, &contacts, "new_vs")?;
        Ok(())
    }

    pub fn enforce_volume_exclusion(
        &mut self,
        old_vs: &[V2d; NVERTS],
        contacts: &[Contact],
    ) -> Result<(), VolExErr> {
        for vi in 0..NVERTS {
            let v = self.poly[vi];
            let old_v = old_vs[vi];
            for contact in contacts {
                self.poly[vi] =
                    fix_point_in_poly(3, old_v, v, &contact.poly)?;
            }
        }

        for vi in 0..NVERTS {
            let wi = circ_ix_plus(vi, NVERTS);
            let v = self.poly[vi];
            let w = self.poly[wi];
            let old_v = old_vs[vi];
            let old_w = old_vs[wi];
            for contact in contacts {
                for other in contact.poly.edges.iter() {
                    let (fixed_v, fixed_w) = fix_edge_intersection(
                        3,
                        (old_v, old_w),
                        (v, w),
                        other,
                    )?;
                    self.poly[vi] = fixed_v;
                    self.poly[wi] = fixed_w;
                }
            }
        }

        Ok(())
    }

    //TODO(BM): automate generation of `num_vars` using proc macro.
    /// Calculate the total number of variables that `CoreState`
    /// holds. That is: the number of variables per vertex, times the
    /// number of all the vertices in a cell.
    pub fn num_vars() -> u32 {
        (NVERTS * 6) as u32
    }

    pub fn flat_sum(&self) -> f64 {
        let mut r: f64 = 0.0;

        for i in 0..(NVERTS) {
            r += self.poly[i].x + self.poly[i].y;
            r += self.rac_acts[i];
            r += self.rac_inacts[i];
            r += self.rho_acts[i];
            r += self.rho_inacts[i];
        }

        r
    }

    pub fn flat_avg(&self) -> f64 {
        self.flat_sum() / (Self::num_vars() as f64)
    }
}

fn violates_volume_exclusion(
    test_v: &V2d,
    test_w: &V2d,
    contacts: &[Contact],
) -> Option<(Poly, V2d, V2d)> {
    for contact in contacts {
        for other in contact.poly.edges.iter() {
            if lsegs_intersect(test_v, test_w, other) {
                return Some((contact.poly, other.p0, other.p1));
            }
        }
    }
    None
}

#[derive(Debug)]
pub enum VolExErr {
    OldEdge,
    OldVert,
    ConfirmViolation,
}

impl From<VolExErr> for String {
    fn from(ve: VolExErr) -> Self {
        format!("{:?}", ve)
    }
}

fn fix_edge_intersection(
    num_iters: usize,
    good_vw: (V2d, V2d),
    new_vw: (V2d, V2d),
    other: &LineSeg2D,
) -> Result<(V2d, V2d), VolExErr> {
    let (orig_v, orig_w) = good_vw;
    if lsegs_intersect(&orig_v, &orig_w, other) {
        return Err(VolExErr::OldEdge);
    }
    let (mut good_v, mut good_w) = good_vw;
    let (mut new_v, mut new_w) = new_vw;
    let mut n = 0;
    while n < num_iters {
        n += 1;
        let test_v = 0.5 * (new_v + good_v);
        let test_w = 0.5 * (new_w + good_w);
        if lsegs_intersect(&test_v, &test_w, other) {
            new_v = test_v;
            new_w = test_w;
        } else {
            good_v = test_v;
            good_w = test_w;
        }
    }
    if lsegs_intersect(&good_v, &good_w, other) {
        Ok((orig_v, orig_w))
    } else {
        Ok((good_v, good_w))
    }
}

fn fix_orig_point(
    orig_v: V2d,
    delta: V2d,
    other: &[V2d; NVERTS],
) -> V2d {
    while is_point_in_poly(&orig_v, None, &other) {
        orig_v = orig_v + delta;
    }
    orig_v
}

fn fix_point_in_poly(
    num_iters: usize,
    mut good_v: V2d,
    mut new_v: V2d,
    other: &Poly,
) -> Result<V2d, VolExErr> {
    let orig_v = good_v;
    if is_point_in_poly(&orig_v, Some(&other.bbox), &other.verts) {
        orig_v = fix_orig_point(
            orig_v,
            (new_v - orig_v).scale(10),
            &other.verts,
        );
    }
    if !is_point_in_poly(&new_v, Some(&other.bbox), &other.verts) {
        return Ok(new_v);
    }
    let mut n = 0;
    while n < num_iters {
        n += 1;
        let test_v = 0.5 * (new_v + good_v);
        if is_point_in_poly(&test_v, Some(&other.bbox), &other.verts)
        {
            new_v = test_v;
        } else {
            good_v = test_v;
        }
    }
    if is_point_in_poly(&good_v, Some(&other.bbox), &other.verts) {
        Ok(orig_v)
    } else {
        Ok(good_v)
    }
}

pub fn confirm_volume_exclusion(
    vs: &[V2d; NVERTS],
    contacts: &[Contact],
    _msg: &str,
) -> Result<(), VolExErr> {
    // use crate::math::v2d::poly_to_string;
    for (vi, v) in vs.iter().enumerate() {
        let wi = circ_ix_plus(vi, NVERTS);
        let w = vs[wi];
        if let Some(_) = violates_volume_exclusion(v, &w, contacts) {
            return Err(VolExErr::ConfirmViolation);
        }
        //     return Err(format!(
        //         "{} violates volume exclusion.\n\
        //             vs = {}, \n\
        //             other poly = {}  \n\
        //             this_vs = {} \n\
        //             a = {}, b = {}",
        //         msg,
        //         v,
        //         &poly_to_string(&p.verts),
        //         &poly_to_string(vs),
        //         a,
        //         b,
        //     ));
        // }
    }
    Ok(())
}
