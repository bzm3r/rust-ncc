#![allow(unused)]
use crate::cell::{CellState, RacRandomState};
use crate::math::{max_f32, min_f32};
use crate::parameters::Parameters;
use crate::world::InteractionState;

type CellDynamicsFn = fn(
    dt: f32,
    state: &CellState,
    rac_random_state: &RacRandomState,
    interaction_state: &InteractionState,
    parameters: &Parameters,
) -> CellState;
const C: [f32; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
// A0s are all zeros
const A1: f32 = 1.0 / 5.0;
const A2: [f32; 2] = [3.0 / 40.0, 9.0 / 40.0];
const A3: [f32; 3] = [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0];
const A4: [f32; 4] = [
    19372.0 / 6561.0,
    -25360.0 / 2187.0,
    64448.0 / 6561.0,
    -212.0 / 729.0,
];
const A5: [f32; 5] = [
    9017.0 / 3168.0,
    -355.0 / 33.0,
    46732.0 / 5247.0,
    49.0 / 176.0,
    -5103.0 / 18656.0,
];
const A6: [f32; 6] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
];
const B: [f32; 7] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
    0.0,
];
const B_HAT: [f32; 7] = [
    5179.0 / 57600.0,
    0.0,
    7571.0 / 16695.0,
    393.0 / 640.0,
    92097.0 / 339200.0,
    187.0 / 2100.0,
    1.0 / 40.0,
];
const INV_QP1: f32 = 1.0 / 5.0; // inverse (max of p and p_hat) + 1, see explanation for equation 4.12 in HNWvol1
const FAC: f32 = 0.8; // safety factor, approximately 0.38^QP1, see explanation for equation 4.12 in HWNvol1
const FAC_MAX: f32 = (5.0 - 1.5) / 2.0; // see explanation for equation 4.12 in HWNvol1

pub struct AuxArgs {
    pub(crate) max_iters: u32,
    pub(crate) atol: f32,
    pub(crate) rtol: f32,
    pub(crate) init_h_factor: Option<f32>,
}

pub struct SolverArgs {
    f: fn(dt: f32, state: &CellState, parameters: &Parameters) -> CellState,
    init_state: CellState,
    t0: f32,
    t1: f32,
}

pub struct Solution {
    pub y: Result<CellState, String>,
    pub num_rejections: u32,
    pub num_iters: u32,
}

pub struct PreMulA {
    // a0s are all zeros
    a1: f32,
    a2: [f32; 2],
    a3: [f32; 3],
    a4: [f32; 4],
    a5: [f32; 5],
    a6: [f32; 6],
}

impl PreMulA {
    fn premul_by(h: f32) -> PreMulA {
        PreMulA {
            a1: h * A1,
            a2: {
                let mut a2 = [0.0_f32; 2];
                (0..2).for_each(|j| a2[j] = h * A2[j]);
                a2
            },
            a3: {
                let mut a3 = [0.0_f32; 3];
                (0..3).for_each(|j| a3[j] = h * A3[j]);
                a3
            },
            a4: {
                let mut a4 = [0.0_f32; 4];
                (0..4).for_each(|j| a4[j] = h * A4[j]);
                a4
            },
            a5: {
                let mut a5 = [0.0_f32; 5];
                (0..5).for_each(|j| a5[j] = h * A5[j]);
                a5
            },
            a6: {
                let mut a6 = [0.0_f32; 6];
                (0..6).for_each(|j| a6[j] = h * A6[j]);
                a6
            },
        }
    }
}

pub struct Ks {
    k0: CellState,
    k1: CellState,
    k2: CellState,
    k3: CellState,
    k4: CellState,
    k5: CellState,
    k6: CellState,
}

impl Ks {
    fn calc(
        premul_ajs: &PreMulA,
        f: CellDynamicsFn,
        dt_primes: [f32; 7],
        init_state: &CellState,
        rand_state: &RacRandomState,
        inter_state: &InteractionState,
        parameters: &Parameters,
    ) -> Ks {
        let PreMulA {
            a1,
            a2,
            a3,
            a4,
            a5,
            a6,
        } = premul_ajs;

        // since dt_primes[0] = 0.0, the function evaluated at that point will return init_state
        let k0 = init_state.clone();

        let k1 = {
            let kp = init_state + a1 * &k0;
            f(dt_primes[1], &kp, rand_state, inter_state, parameters)
        };

        let k2 = {
            let kp = init_state + a2[0] * &k0 + a2[1] * &k1;
            f(dt_primes[2], &kp, rand_state, inter_state, parameters)
        };

        let k3 = {
            let kp = init_state + a3[0] * &k0 + a3[1] * &k1 + a3[2] * &k2;
            f(dt_primes[3], &kp, rand_state, inter_state, parameters)
        };

        let k4 = {
            let kp = init_state + a4[0] * &k0 + a4[1] * &k1 + a4[2] * &k2 + a4[3] * &k3;
            f(dt_primes[4], &kp, rand_state, inter_state, parameters)
        };

        let k5 = {
            let kp =
                init_state + a5[0] * &k0 + a5[1] * &k1 + a5[2] * &k2 + a5[3] * &k3 + a5[4] * &k4;
            f(dt_primes[5], &kp, rand_state, inter_state, parameters)
        };

        let k6 = {
            let kp = init_state
                + a6[0] * &k0
                + a6[1] * &k1
                + a6[2] * &k2
                + a6[3] * &k3
                + a6[4] * &k4
                + a6[5] * &k5;
            f(dt_primes[6], &kp, rand_state, inter_state, parameters)
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

pub fn rkdp5(
    mut dt: f32,
    f: CellDynamicsFn,
    mut init_state: &CellState,
    rand_state: &RacRandomState,
    inter_state: &InteractionState,
    parameters: &Parameters,
    mut aux_args: AuxArgs,
) -> Solution {
    let mut init_state = init_state.clone();

    let AuxArgs {
        max_iters,
        atol,
        rtol,
        init_h_factor,
    } = aux_args;

    let mut h = if let Some(h_factor) = init_h_factor {
        h_factor * dt
    } else {
        0.1 * dt
    };

    let mut num_iters = 0_u32;
    let mut fac_max = FAC_MAX;
    let mut last_iter = false;
    let mut num_rejections: u32 = 0;

    while num_iters < max_iters {
        // premultiply Ajs (j in 0..6) by h
        let premul_as = PreMulA::premul_by(h);

        let dt_primes = {
            let mut dt_primes = [0.0_f32; 7];
            (0..7).for_each(|j| dt_primes[j] = h * C[j]);
            dt_primes
        };

        let Ks {
            k0,
            k1,
            k2,
            k3,
            k4,
            k5,
            k6,
        } = Ks::calc(
            &premul_as,
            f,
            dt_primes,
            &init_state,
            rand_state,
            inter_state,
            parameters,
        );

        let PreMulA { a6, .. } = premul_as;
        // because B is the same as A6, except with last value 0, we can again use the premultiplied values of A6
        let y1 = init_state
            + a6[0] * k0
            + a6[1] * k1
            + a6[2] * k2
            + a6[3] * k3
            + a6[4] * k4
            + a6[5] * k5;

        if last_iter {
            assert!((h - dt).abs() < f32::EPSILON);
            return Solution {
                y: Ok(y1),
                num_rejections,
                num_iters,
            };
        }

        // on the other hand, we don't have pre-multiplied values for B_HAT
        let y1_hat = init_state
            + h * (B_HAT[0] * k0
                + B_HAT[1] * k1
                + B_HAT[2] * k2
                + B_HAT[3] * k3
                + B_HAT[4] * k4
                + B_HAT[5] * k5
                + B_HAT[6] * k6);

        // Equations 4.10, 4.11, Hairer,Wanner&Norsett Solving ODEs Vol. 1
        let sc = init_state
            .abs()
            .max(&y1.abs())
            .scalar_mul(rtol)
            .scalar_add(atol);
        let delta = y1 - y1_hat;
        let error = (delta.powi(2) / sc).average().sqrt();
        let mut h_new = h * min_f32(fac_max, FAC * (1.0 / error).powf(INV_QP1));

        // see explanation for equation 4.13 in HNWvol1
        if error <= 1.0 {
            fac_max = FAC_MAX;
            init_state = y1;
            if h_new > dt - h {
                h_new = dt - h;
                last_iter = true;
            };
            dt -= h;
            h = h_new;
        } else {
            num_rejections += 1;
            h = h_new;
            fac_max = 1.0;
        }

        num_iters += 1;
    }

    Solution {
        num_rejections,
        num_iters,
        y: Err(format!("Too many iterations!")),
    }
}
