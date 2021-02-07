import experiment_templates as ets
import hardio as fw
import numpy as np

TIME_IN_HOURS = 0.01
T = 2.0
L = 1e-6
F = 1e-9
ETA = 0.1
L3D = 10e-6
K_MEM_ON = 0.02
K_MEM_OFF = 0.15
KGTP = 1e-4
COA = 0
CIL = 0
INIT_RAC = np.array([0.3 / 4.0] * 4 + [0.0] * (16 - 4))

NUM_CELLS = 2
BOX_WIDTH = 2
BOX_HEIGHT = 1
NUM_INT_STEPS = 10
TIME_IN_SECS = TIME_IN_HOURS * 3600.0
NUM_TSTEPS = int(TIME_IN_SECS / T)
NVERTS = 16
RAND_SCHEME = None


def u(x):
    if x > 0:
        return 1.0
    else:
        return 0.0


params = dict([
    ("16", NVERTS),
    ("num_tsteps", NUM_TSTEPS),
    ("num_cells", NUM_CELLS),
    ("num_int_steps", NUM_INT_STEPS),
    ("t", T),
    ("l", L),
    ("f", F),
    ("eta", ETA),
    ("l3d", L3D),
    ("k_mem_on", K_MEM_ON),
    ("k_mem_off", K_MEM_OFF),
    ("kgtp", KGTP),
    ("kdgtp", KGTP),
    ("cil_mag", CIL),
    ("coa_mag", COA),
    ("coa_los_penalty", 2.0 * u(COA)),
    ("coa_range", 220e-6 * u(COA)),
    ("close_zero_at", 1.5e-6),
    ("close_one_at", 0.5e-6),
    ("randomization", False),
    ("cell_r", 20e-6),
    ("tot_rac", 2.5e6),
    ("tot_rho", 1e6),
    ("init_act_rgtp", 0.1),
    ("init_inact_rgtp", 0.1),
    ("diffusion_rgtp", 0.1e-12),
    ("kgdi", 1),
    ("kdgdi", 1),
    ("halfmax_vertex_rgtp_act", 0.4),
    ("threshold_rho_activity", 0.4),
    ("stiffness_cyto", 1e-5),
    ("force_rho", 0.2),
    ("rand_avg_t", 40.0),
    ("rand_std_t", 0.1 * 40.0),
    ("rand_mag", 10.0),
    ("kgtp_rac", 24.0),
    ("kgtp_rho", 28.0),
    ("kdgtp_rac", 8.0),
    ("kdgtp_rho", 60.0),
    ("kgtp_rac_autoact", 500.0),
    ("kgtp_rho_autoact", 390.0),
    ("kdgtp_rac_on_rho", 400.0),
    ("kdgtp_rho_on_rac", 4000.0),
    ("halfmax_tension_inhib", 0.1),
    ("tension_inhib", 40.0),
    ("const_protrusive", 3000.0),
    ("stiffness_edge", 8000.0),
    ("lm_h", 200e-9),
    ("lm_ss", 10e3),
    ("halfmax_rgtp_max_f_frac", 0.3),
    ("num_rand_vs", 16 * 0.25),
    ("total_rgtp", 1.0),
    ("init_rac", INIT_RAC),
    ("init_rho", np.roll(INIT_RAC, 8)),
    ("vertex_eta", 2.9),
])

# coa_dict = {
#     49: 8.0,
#     36: 9.0,
#     25: 12.0,
#     16: 14.0,
#     9: 16.0,
#     4: 24.0,
#     2: 24.0,
#     1: 24.0}

if __name__ == "__main__":
    uniform_initial_polarization = False

    ets.rust_comparison_test(
        params,
        box_width=BOX_WIDTH,
        box_height=BOX_HEIGHT,
    )
