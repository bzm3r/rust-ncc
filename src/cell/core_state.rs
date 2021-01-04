use crate::cell::chemistry::{
    calc_conc_rgtps, calc_kdgtps_rac, calc_kdgtps_rho,
    calc_kgtps_rac, calc_kgtps_rho, calc_net_fluxes, RacRandState,
    RgtpDistribution,
};
use crate::cell::mechanics::{
    calc_cyto_forces, calc_edge_forces, calc_edge_vecs,
    calc_rgtp_forces,
};
use crate::interactions::{CellInteractions, RgtpActivityDiff};
use crate::math::v2d::V2D;
use crate::math::{hill_function3, max_f32, min_f32};
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::circ_ix_minus;
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

/// `CoreState` contains all the variables that are simulated within
/// geometric updates. That is, they are modelled using ODEs which
/// are then integrated using either the Euler method or
/// Runge-Kutta Dormand-Prince 5 (Matlab's `ode45`).
#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize)]
pub struct CoreState {
    /// Polygon representing cell shape.
    pub poly: [V2D; NVERTS],
    /// Fraction of Rac1 active at each vertex.
    rac_acts: [f32; NVERTS],
    /// Fraction of Rac1 inactive at each vertex.
    rac_inacts: [f32; NVERTS],
    /// Fraction of RhoA active at each vertex.
    rho_acts: [f32; NVERTS],
    /// Fraction of RhoA inactive at each vertex.
    rho_inacts: [f32; NVERTS],
}

impl Add for CoreState {
    type Output = CoreState;

    fn add(self, rhs: CoreState) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = self.poly[i] + rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] + rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] + rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] + rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] + rhs.rho_inacts[i]
        }

        Self::Output {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Sub for CoreState {
    type Output = CoreState;

    fn sub(self, rhs: CoreState) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = self.poly[i] - rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] - rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] - rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] - rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] - rhs.rho_inacts[i]
        }

        Self::Output {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Div for CoreState {
    type Output = CoreState;

    fn div(self, rhs: CoreState) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = self.poly[i] / rhs.poly[i];
            rac_acts[i] = self.rac_acts[i] / rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] / rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] / rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] / rhs.rho_inacts[i]
        }

        Self::Output {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Mul<CoreState> for f32 {
    type Output = CoreState;

    fn mul(self, rhs: CoreState) -> CoreState {
        let mut poly = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            poly[i] = self * rhs.poly[i];
            rac_acts[i] = self * rhs.rac_acts[i];
            rac_inacts[i] = self * rhs.rac_inacts[i];
            rho_acts[i] = self * rhs.rho_acts[i];
            rho_inacts[i] = self * rhs.rho_inacts[i]
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

/// Records the mechanical state of a cell.
#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize)]
pub struct MechState {
    /// Strain each edge is under, where resting edge length is
    /// defined in the cell's parameters.
    pub edge_strains: [f32; NVERTS],
    /// Forces on each vertex due to Rho GTPase activity.
    pub rgtp_forces: [V2D; NVERTS],
    /// Forces on each vertex due to cytoplasmic pressure.
    pub cyto_forces: [V2D; NVERTS],
    /// Forces on each vertex due to edge-edge (elastic) forces.
    pub edge_forces: [V2D; NVERTS],
    /// Average of the strain in edges which are under tension (i.e.
    /// they are longer than their initial resting edge length.
    pub avg_tens_strain: f32,
    /// Sum of all forces that are acting on a vertex, except for
    /// adhesion, which comes from interaction information.
    pub sum_fs: [V2D; NVERTS],
}

/// Calculates the various rates necessary to define the ODEs
/// simulating biochemistry. These are (`X` is either `rac` or `rho`):
///     * `kdgtp_X`: Rho GTPase inactivation rates
///     * `kgtp_X`: Rho GTPase activation rates
///     * `X_act_net_fluxes`: diffusion fluxes between vertices of
/// active form of Rho GTPases
///     * `X_inact_net_fluxes`: diffusion fluxes between vertices of
/// inactive form of Rho GTPases
///     * `X_cyto`: fraction of Rho GTPase in the cytoplasm
///     * `x_tens`: "tension" factor that affects Rac1 activation
/// rate, calculated based on average tensile strain in cell (i.e.
/// how stretched the cell is).
#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize)]
pub struct ChemState {
    pub kdgtps_rac: [f32; NVERTS],
    pub kgtps_rac: [f32; NVERTS],
    pub rac_act_net_fluxes: [f32; NVERTS],
    pub rac_inact_net_fluxes: [f32; NVERTS],
    pub kdgtps_rho: [f32; NVERTS],
    pub kgtps_rho: [f32; NVERTS],
    pub rac_cyto: f32,
    pub rho_cyto: f32,
    pub rho_act_net_fluxes: [f32; NVERTS],
    pub rho_inact_net_fluxes: [f32; NVERTS],
    pub x_tens: f32,
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize)]
pub struct GeomState {
    /// Unit edge vectors which point from position of vertex `vi`
    /// to position of vertex `vi + 1`, where `vi + 1` is calculated
    /// modulo `NVERTS`. Note that an edge is defined by its "lower"
    /// index, modulo `NVERTS`. That is, the edge `(0, 1)`, is
    /// different from the edge `(1, 2)`, and `(0, 1)` is also
    /// different from `(15, 0)` (assuming that `NVERTS == 16` in
    /// this example).
    pub unit_edge_vecs: [V2D; NVERTS],
    /// Length of edges. Each edge is defined by its smallest vertex.
    pub edge_lens: [f32; NVERTS],
    /// Inward pointing unit vectors at each vertex. These are
    /// calculated so that they bisect the angle between the two
    /// edges which meet at a vertex.
    pub unit_inward_vecs: [V2D; NVERTS],
}

// #[derive(Copy, Clone, Debug, Default, Deserialize, Serialize)]
// pub struct DepStates {
//     pub geom_state: GeomState,
//     pub chem_state: ChemState,
//     pub mech_state: MechState,
// }

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

impl Display for CoreState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // println!("----------");
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
        fmt_var_arr(f, "uivs", &self.unit_inward_vecs)?;
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
        fmt_var_arr(f, "tot_forces", &self.sum_fs)
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

// impl Display for DepStates {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         writeln!(f, "{}", &self.geom_state)?;
//         writeln!(f, "{}", &self.mech_state)?;
//         writeln!(f, "{}", &self.chem_state)
//     }
// }

impl CoreState {
    //TODO(BM): automate generation of `num_vars` using proc macro.
    /// Calculate the total number of variables that `CoreState`
    /// holds. That is: the number of variables per vertex, times the
    /// number of all the vertices in a cell.
    pub fn num_vars() -> usize {
        (NVERTS * 6) as usize
    }

    pub fn calc_geom_state(&self) -> GeomState {
        // Calculate edge vectors of a polygon.
        let evs = calc_edge_vecs(&self.poly);
        // Calculate magnitude of each edge vec, to get its length.
        let mut edge_lens = [0.0_f32; NVERTS];
        (0..NVERTS).for_each(|i| edge_lens[i] = (&evs[i]).mag());
        // Divide each edge vector by its magnitude to get the
        // corresponding unit vector.
        let mut uevs = [V2D::default(); NVERTS];
        (0..NVERTS).for_each(|i| uevs[i] = (&evs[i]).unitize());
        // Given two unit edge vectors, find the vector which points
        // into the polygon and bisects the angle
        let mut uivs = [V2D::default(); NVERTS];
        (0..NVERTS).for_each(|i| {
            let im1 = circ_ix_minus(i as usize, NVERTS);
            let tangent = (uevs[i] + uevs[im1]).unitize();
            uivs[i] = tangent.normal();
        });

        GeomState {
            unit_edge_vecs: uevs,
            edge_lens,
            unit_inward_vecs: uivs,
        }
    }

    pub fn calc_mech_state(
        &self,
        geom_state: &GeomState,
        parameters: &Parameters,
    ) -> MechState {
        let GeomState {
            unit_edge_vecs: uevs,
            edge_lens,
            unit_inward_vecs: uivs,
        } = geom_state;
        let rgtp_forces = calc_rgtp_forces(
            &self.rac_acts,
            &self.rho_acts,
            uivs,
            parameters.halfmax_vertex_rgtp_act,
            parameters.const_protrusive,
            parameters.const_retractive,
        );
        let cyto_forces = calc_cyto_forces(
            &self.poly,
            &uivs,
            parameters.rest_area,
            parameters.stiffness_ctyo,
        );
        // Calculate strain in each edge.
        let mut edge_strains = [0.0_f32; NVERTS];
        (0..NVERTS).for_each(|i| {
            edge_strains[i] =
                (edge_lens[i] / parameters.rest_edge_len) - 1.0
        });
        let edge_forces = calc_edge_forces(
            &edge_strains,
            uevs,
            parameters.stiffness_edge,
        );
        // If strain is positive (tensile), then consider it in this
        // averaging, otherwise don't. This is because I'm assuming
        // that compression does not have an effect on Rac1 activity.
        // Only tension is considered to have an effect (see refs. in
        // SI.
        //TODO(BM): what is the latest on this front? Recent paper
        // (Elife?) which suggests not true for migrating cells.
        let avg_tens_strain = edge_strains
            .iter()
            .map(|&es| if es < 0.0 { 0.0 } else { es })
            .sum::<f32>()
            / NVERTS as f32;
        // Sum of all the non-adhesive forces acting on the cell.
        let mut sum_fs = [V2D::default(); NVERTS];
        (0..NVERTS).for_each(|i| {
            sum_fs[i] =
                rgtp_forces[i] + cyto_forces[i] + edge_forces[i]
                    - edge_forces[circ_ix_minus(i as usize, NVERTS)];
        });
        MechState {
            edge_strains,
            rgtp_forces,
            cyto_forces,
            edge_forces,
            avg_tens_strain,
            sum_fs,
        }
    }

    pub fn calc_chem_state(
        &self,
        geom_state: &GeomState,
        mech_state: &MechState,
        rac_rand_state: &RacRandState,
        inter_state: &CellInteractions,
        parameters: &Parameters,
    ) -> ChemState {
        let GeomState { edge_lens, .. } = geom_state;
        // Need to calculate average length of edges meeting at
        // a vertex in order to roughly approximate diffusion related
        // flux of Rho GTPases from neighbouring vertices. Provides
        // an approximation for the length of the membrane abstracted
        // by the edges that meet at that vertex.
        let mut avg_edge_lens: [f32; NVERTS] = [0.0_f32; NVERTS];
        (0..NVERTS).for_each(|i| {
            let im1 = circ_ix_minus(i as usize, NVERTS);
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
            &inter_state.x_coas,
            &inter_state.x_chem_attrs,
            &inter_state.x_cals,
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
            &inter_state.x_cils,
            x_tens,
            parameters.kdgtp_rac,
            parameters.kdgtp_rho_on_rac,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let kgtps_rho = calc_kgtps_rho(
            &self.rho_acts,
            &conc_rho_acts,
            &inter_state.x_cils,
            parameters.kgtp_rho,
            parameters.halfmax_vertex_rgtp_conc,
            parameters.kgtp_rho_auto,
        );
        let kdgtps_rho = calc_kdgtps_rho(
            &self.rho_acts,
            &conc_rac_acts,
            &inter_state.x_cals,
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

        let rac_cyto = parameters.total_rgtp
            - self.rac_acts.iter().sum::<f32>()
            - self.rac_inacts.iter().sum::<f32>();
        let rho_cyto = parameters.total_rgtp
            - self.rho_acts.iter().sum::<f32>()
            - self.rho_inacts.iter().sum::<f32>();
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
    /// vertex motion and biochemistry. In particular, calculate
    /// `(delta(state)/delta(t))`. Note that `delta(state)` should have
    /// "units" of `[CoreState]` (the units of the result of
    /// addition/subtraction of two quanties with units X, is X), so
    /// this function should be returning a quantity with units
    /// `[CoreState]/[Time]`, but since time is normalized, this is
    /// the same as having units of `[CoreState]`.
    pub fn dynamics_f(
        state: &CoreState,
        rac_rand_state: &RacRandState,
        inter_state: &CellInteractions,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> CoreState {
        let geom_state = Self::calc_geom_state(state);
        let mech_state =
            Self::calc_mech_state(state, &geom_state, parameters);
        let chem_state = Self::calc_chem_state(
            state,
            &geom_state,
            &mech_state,
            rac_rand_state,
            inter_state,
            parameters,
        );
        let mut delta = CoreState::default();
        for i in 0..NVERTS {
            // rate of rac deactivation * current fraction of rac active
            let inactivated_rac =
                chem_state.kdgtps_rac[i] * state.rac_acts[i];
            let activated_rac =
                chem_state.kgtps_rac[i] * state.rac_inacts[i];
            let delta_rac_activated = activated_rac - inactivated_rac;
            let rac_cyto_exchange = {
                let rac_mem_on =
                    parameters.k_mem_on_vertex * chem_state.rac_cyto;
                let rac_mem_off =
                    parameters.k_mem_off * state.rac_inacts[i];
                rac_mem_on - rac_mem_off
            };
            let vertex_rac_act_flux =
                chem_state.rac_act_net_fluxes[i];
            let vertex_rac_inact_flux =
                chem_state.rac_inact_net_fluxes[i];
            delta.rac_acts[i] =
                delta_rac_activated + vertex_rac_act_flux;
            delta.rac_inacts[i] = rac_cyto_exchange
                + vertex_rac_inact_flux
                - delta_rac_activated;

            let inactivated_rho =
                chem_state.kdgtps_rho[i] * state.rho_acts[i];
            let activated_rho =
                chem_state.kgtps_rho[i] * state.rho_inacts[i];
            let delta_rho_activated = activated_rho - inactivated_rho;
            let rho_cyto_exchange = {
                let rho_mem_on =
                    parameters.k_mem_on_vertex * chem_state.rho_cyto;
                let rho_mem_off =
                    parameters.k_mem_off * state.rho_inacts[i];
                rho_mem_on - rho_mem_off
            };
            let vertex_rho_act_flux =
                chem_state.rho_act_net_fluxes[i];
            let vertex_rho_inact_flux =
                chem_state.rho_inact_net_fluxes[i];
            delta.rho_acts[i] =
                delta_rho_activated + vertex_rho_act_flux;
            delta.rho_inacts[i] = rho_cyto_exchange
                + vertex_rho_inact_flux
                - delta_rho_activated;
            delta.poly[i] = (1.0 / world_parameters.vertex_eta)
                * (mech_state.sum_fs[i] + inter_state.x_adhs[i]);
        }
        delta
    }

    pub fn new(
        vertex_coords: [V2D; NVERTS],
        init_rac: RgtpDistribution,
        init_rho: RgtpDistribution,
    ) -> CoreState {
        // x_cils: [f32; NVERTS], x_coas: [f32; NVERTS], x_chemoas: [f32; NVERTS], x_rands: [f32; NVERTS], x_bdrys: [f32; NVERTS];
        CoreState {
            poly: vertex_coords,
            rac_acts: init_rac.active,
            rac_inacts: init_rac.inactive,
            rho_acts: init_rho.active,
            rho_inacts: init_rho.inactive,
        }
    }

    pub fn scalar_mul(&self, s: f32) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = s * self.poly[i];
            rac_acts[i] = self.rac_acts[i] * s;
            rac_inacts[i] = self.rac_inacts[i] * s;
            rho_acts[i] = self.rho_acts[i] * s;
            rho_inacts[i] = self.rho_inacts[i] * s;
        }

        CoreState {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn scalar_add(&self, s: f32) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = s + self.poly[i];
            rac_acts[i] = self.rac_acts[i] + s;
            rac_inacts[i] = self.rac_inacts[i] + s;
            rho_acts[i] = self.rho_acts[i] + s;
            rho_inacts[i] = self.rho_inacts[i] + s;
        }

        CoreState {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn abs(&self) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].abs();
            rac_acts[i] = self.rac_acts[i].abs();
            rac_inacts[i] = self.rac_inacts[i].abs();
            rho_acts[i] = self.rho_acts[i].abs();
            rho_inacts[i] = self.rho_inacts[i].abs();
        }

        CoreState {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn powi(&self, x: i32) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].powi(x);
            rac_acts[i] = self.rac_acts[i].powi(x);
            rac_inacts[i] = self.rac_inacts[i].powi(x);
            rho_acts[i] = self.rho_acts[i].powi(x);
            rho_inacts[i] = self.rho_inacts[i].powi(x);
        }

        CoreState {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn max(&self, other: &CoreState) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].max(&other.poly[i]);
            rac_acts[i] =
                max_f32(self.rac_acts[i], other.rac_acts[i]);
            rac_inacts[i] =
                max_f32(self.rac_inacts[i], other.rac_inacts[i]);
            rho_acts[i] =
                max_f32(self.rho_acts[i], other.rho_acts[i]);
            rho_inacts[i] =
                max_f32(self.rho_inacts[i], other.rho_inacts[i]);
        }

        CoreState {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn min(&self, other: &CoreState) -> CoreState {
        let mut vertex_coords = [V2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].min(&other.poly[i]);
            rac_acts[i] =
                min_f32(self.rac_acts[i], other.rac_acts[i]);
            rac_inacts[i] =
                min_f32(self.rac_inacts[i], other.rac_inacts[i]);
            rho_acts[i] =
                min_f32(self.rho_acts[i], other.rho_acts[i]);
            rho_inacts[i] =
                min_f32(self.rho_inacts[i], other.rho_inacts[i]);
        }

        CoreState {
            poly: vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn sum(&self) -> f32 {
        let mut r: f32 = 0.0;

        for i in 0..(NVERTS) {
            r += self.poly[i].x + self.poly[i].y;
            r += self.rac_acts[i];
            r += self.rac_inacts[i];
            r += self.rho_acts[i];
            r += self.rho_inacts[i];
        }

        r
    }

    pub fn average(&self) -> f32 {
        self.sum() / (Self::num_vars() as f32)
    }

    /// Calculate which Rho GTPase has dominates in terms of effect
    /// at this vertex.
    pub fn calc_crl_rgtp_state(
        &self,
        parameters: &Parameters,
    ) -> [RgtpActivityDiff; NVERTS] {
        let mut r = [0.0; NVERTS];
        self.rac_acts
            .iter()
            .zip(self.rho_acts.iter())
            .enumerate()
            .for_each(|(ix, (&rac, &rho))| {
                r[ix] = hill_function3(
                    parameters.halfmax_vertex_rgtp_act,
                    rac,
                ) - hill_function3(
                    parameters.halfmax_vertex_rgtp_act,
                    rho,
                );
            });
        r
    }

    #[cfg(feature = "debug_mode")]
    pub fn validate(
        &self,
        loc_str: &str,
        parameters: &Parameters,
    ) -> Result<(), String> {
        if self.rac_acts.iter().any(|&r| r < 0.0_f32) {
            return Err(format!(
                "{}: neg rac_acts: {:?}",
                loc_str, self.rac_acts
            ));
        }
        if self.rac_inacts.iter().any(|&r| r < 0.0_f32) {
            return Err(format!(
                "{}: neg rac_inacts: {:?}",
                loc_str, self.rac_inacts
            ));
        }
        if self.rho_acts.iter().any(|&r| r < 0.0_f32) {
            return Err(format!(
                "{}: neg rho_acts: {:?}",
                loc_str, self.rho_acts
            ));
        }
        if self.rho_inacts.iter().any(|&r| r < 0.0_f32) {
            return Err(format!(
                "{}: neg rho_inacts: {:?}",
                loc_str, self.rho_inacts
            ));
        }
        let sum_rac_mem = self.rac_inacts.iter().sum::<f32>()
            + self.rac_acts.iter().sum::<f32>();
        if sum_rac_mem > parameters.total_rgtp || sum_rac_mem < 0.0 {
            return Err(format!(
                "{}: problem in sum of rac_mem: {}",
                loc_str, sum_rac_mem
            ));
        }
        let sum_rho_mem = self.rho_inacts.iter().sum::<f32>()
            + self.rho_acts.iter().sum::<f32>();
        if sum_rho_mem > parameters.total_rgtp || sum_rho_mem < 0.0 {
            return Err(format!(
                "{}: problem in sum of rho_mem: {}",
                loc_str, sum_rho_mem
            ));
        }
        Ok(())
        // println!("{}: successfully validated", loc_str)
    }
}
