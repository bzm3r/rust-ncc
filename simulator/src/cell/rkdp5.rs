use crate::cell::states::DCoreDt;
use crate::cell::{chemistry::RacRandState, states::Core};
use crate::interactions::{ContactData, Interactions};
use crate::math::min_f64;
use crate::parameters::{Parameters, WorldParameters};
use crate::world::RkOpts;

type CellDynamicsFn = fn(
    state: &Core,
    rac_random_state: &RacRandState,
    interactions: &Interactions,
    world_parameters: &WorldParameters,
    parameters: &Parameters,
) -> DCoreDt;

// A0s are all zeros
const A1: f64 = 1.0 / 5.0;
const A2: [f64; 2] = [3.0 / 40.0, 9.0 / 40.0];
const A3: [f64; 3] = [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0];
const A4: [f64; 4] = [
    19372.0 / 6561.0,
    -25360.0 / 2187.0,
    64448.0 / 6561.0,
    -212.0 / 729.0,
];
const A5: [f64; 5] = [
    9017.0 / 3168.0,
    -355.0 / 33.0,
    46732.0 / 5247.0,
    49.0 / 176.0,
    -5103.0 / 18656.0,
];
const A6: [f64; 6] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
];
const B: [f64; 7] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
    0.0,
];
const B_HAT: [f64; 7] = [
    5179.0 / 57600.0,
    0.0,
    7571.0 / 16695.0,
    393.0 / 640.0,
    92097.0 / 339200.0,
    187.0 / 2100.0,
    1.0 / 40.0,
];
const INV_QP1: f64 = 1.0 / 5.0; // inverse (max of p and p_hat) + 1, see explanation for equation 4.12 in HNW vol1
const FAC: f64 = 0.8; // safety factor, approximately 0.38^QP1, see explanation for equation 4.12 in HWN vol1
const FAC_MAX: f64 = (5.0 - 1.5) / 2.0; // see explanation for equation 4.12 in HWN vol1

pub struct Solution {
    pub state: Result<Core, String>,
    pub num_rejections: usize,
    pub num_iters: usize,
}

pub struct Ks {
    k0: DCoreDt,
    k1: DCoreDt,
    k2: DCoreDt,
    k3: DCoreDt,
    k4: DCoreDt,
    k5: DCoreDt,
    k6: DCoreDt,
}

impl Ks {
    fn calc(
        f: CellDynamicsFn,
        h: f64,
        init_state: Core,
        rand_state: &RacRandState,
        inter_state: &Interactions,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> Ks {
        let k0 = f(
            &init_state,
            rand_state,
            inter_state,
            world_parameters,
            parameters,
        );

        let k1 = {
            let kp = init_state + h * k0.time_step(A1);
            f(
                &kp,
                rand_state,
                inter_state,
                world_parameters,
                parameters,
            )
        };

        let k2 = {
            let kp = init_state
                + h * (k0.time_step(A2[0]) + k1.time_step(A2[1]));
            f(
                &kp,
                rand_state,
                inter_state,
                world_parameters,
                parameters,
            )
        };

        let k3 = {
            let kp = init_state
                + h * (k0.time_step(A3[0])
                    + k1.time_step(A3[1])
                    + k2.time_step(A3[2]));
            f(
                &kp,
                rand_state,
                inter_state,
                world_parameters,
                parameters,
            )
        };

        let k4 = {
            let kp = init_state
                + h * (k0.time_step(A4[0])
                    + k1.time_step(A4[1])
                    + k2.time_step(A4[2])
                    + k3.time_step(A4[3]));
            f(
                &kp,
                rand_state,
                inter_state,
                world_parameters,
                parameters,
            )
        };

        let k5 = {
            let kp = init_state
                + h * (k0.time_step(A5[0])
                    + k1.time_step(A5[1])
                    + k2.time_step(A5[2])
                    + k3.time_step(A5[3])
                    + k4.time_step(A5[4]));
            f(
                &kp,
                rand_state,
                inter_state,
                world_parameters,
                parameters,
            )
        };

        let k6 = {
            let kp = init_state
                + h * (k0.time_step(A6[0])
                    + k1.time_step(A6[1])
                    + k2.time_step(A6[2])
                    + k3.time_step(A6[3])
                    + k4.time_step(A6[4])
                    + k5.time_step(A6[5]));
            f(
                &kp,
                rand_state,
                inter_state,
                world_parameters,
                parameters,
            )
        };

        Ks {
            k0,
            k1,
            k2,
            k3,
            k4,
            k5,
            k6,
        }
    }
}

pub fn integrator(
    mut dt: f64,
    f: CellDynamicsFn,
    mut init_state: Core,
    rand_state: &RacRandState,
    interactions: &Interactions,
    world_parameters: &WorldParameters,
    parameters: &Parameters,
    contact_data: &[ContactData],
    int_opts: RkOpts,
) -> Solution {
    let RkOpts {
        max_iters,
        atol,
        rtol,
        init_h_scale,
        ..
    } = int_opts;

    let mut h = init_h_scale * dt;

    let mut num_iters: usize = 0;
    let mut fac_max = FAC_MAX;
    let mut last_iter = false;
    let mut num_rejections: usize = 0;

    while num_iters < max_iters {
        let Ks {
            k0,
            k1,
            k2,
            k3,
            k4,
            k5,
            k6,
        } = Ks::calc(
            f,
            h,
            init_state,
            rand_state,
            interactions,
            world_parameters,
            parameters,
        );

        let mut next_state = init_state
            + B[0] * k0.time_step(h)
            + B[1] * k1.time_step(h)
            + B[2] * k2.time_step(h)
            + B[3] * k3.time_step(h)
            + B[4] * k4.time_step(h)
            + B[5] * k5.time_step(h)
            + B[6] * k6.time_step(h);

        if last_iter {
            assert!((h - dt).abs() < f64::EPSILON);
            next_state
                .strict_enforce_volume_exclusion(
                    &init_state.poly,
                    &contact_data,
                )
                .map_or_else(|e| panic!(e), |_| {});
            return Solution {
                state: Ok(next_state),
                num_rejections,
                num_iters,
            };
        }

        let next_state_hat = init_state
            + h * (k0.time_step(B_HAT[0])
                + k1.time_step(B_HAT[1])
                + k2.time_step(B_HAT[2])
                + k3.time_step(B_HAT[3])
                + k4.time_step(B_HAT[4])
                + k5.time_step(B_HAT[5])
                + k6.time_step(B_HAT[6]));

        // Equations 4.10, 4.11, Hairer,Wanner&Norsett Solving ODEs Vol. 1
        let sc =
            rtol * init_state.abs().max(&next_state.abs()) + atol;
        let error = ((next_state - next_state_hat).square() / sc)
            .flat_avg()
            .sqrt();
        let mut h_new =
            h * min_f64(fac_max, FAC * (1.0 / error).powf(INV_QP1));

        // see explanation for equation 4.13 in HNW vol1
        if error <= 1.0 {
            fac_max = FAC_MAX;
            next_state
                .strict_enforce_volume_exclusion(
                    &init_state.poly,
                    &contact_data,
                )
                .map_or_else(|e| panic!(e), |_| {});
            init_state = next_state;
            if h + h_new > dt {
                h_new = dt - h;
                last_iter = true;
            };
            dt -= h;
            h = h_new;
        } else {
            fac_max = 1.0;
            num_rejections += 1;
            h = h_new;
        }

        num_iters += 1;
    }

    Solution {
        state: Err("Too many iterations!".to_string()),
        num_rejections,
        num_iters,
    }
}
