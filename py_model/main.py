import experiment_templates as ets
import hardio as fw
import numpy as np
import argparse
import os

print("Simulating with Python model.")
parser = argparse.ArgumentParser()
parser.add_argument("--name", help="name", type=str)
parser.add_argument("--final_t",
                    help="simulation run time, in seconds",
                    type=float)
# parser.add_argument("--num_int_steps",
#                     help="number of intermediate integration "
#                          "steps in the forward euler timestepper",
#                     type=float)
parser.add_argument("--snap_period", help="period between saved snapshots",
                    type=float)
parser.add_argument("--num_cells",
                    help="number of cells in simulation",
                    type=int)
# parser.add_argument("--arrangement",
#                     help="a string specifying the name of the arrangement "
#                          "according to which the simulation is initially set",
#                     type=int)
# parser.add_argument("--box_height",
#                     help="height of box in which cells are initially placed",
#                     type=int)
parser.add_argument("--cil",
                    help="cil mag",
                    type=int)
parser.add_argument("--coa",
                    help="coa mag",
                    type=int)
parser.add_argument("--out_dir", help="python output dir", type=str)

args = parser.parse_args()

T = 2.0
NAME = args.name
FINAL_T = args.final_t
SNAP_PERIOD = args.snap_period
NUM_INT_STEPS = 10
NUM_CELLS = args.num_cells
BOX_WIDTH = 2
BOX_HEIGHT = 1
COA = args.coa
CIL = args.cil
OUT_DIR = os.path.normpath(args.out_dir)
if not os.path.exists(OUT_DIR):
    print("creating non-existent out dir: {}".format(OUT_DIR))
    os.mkdir(OUT_DIR)
else:
    print("out dir: {}".format(OUT_DIR))

L = 1e-6
F = 1e-9
ETA = 0.1
L3D = 10e-6
K_MEM_ON = 0.02
K_MEM_OFF = 0.15
KGTP = 1e-4
INIT_RAC = np.array([0.3 / 4.0] * 4 + [0.0] * (16 - 4))
NVERTS = 16
RAND_SCHEME = None

LOG_LEVEL = fw.NO_LOG


def u(x):
    if x > 0:
        return 1.0
    else:
        return 0.0


params = dict([
    ("16", NVERTS),
    ("final_t", FINAL_T),
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
    ("cil_mag", CIL),
    ("coa_mag", COA),
    ("coa_los_penalty", 2.0 * u(COA)),
    ("coa_halfmax_dist", 0.5 * 220e-6 * u(COA)),
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
    ("halfmax_vertex_rgtp", 0.4),
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
    ("snap_period", SNAP_PERIOD),
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
        OUT_DIR,
        NAME,
        99,
        box_width=BOX_WIDTH,
        box_height=BOX_HEIGHT,
    )
