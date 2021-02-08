from retrieve import *
import numpy as np
import matplotlib.pyplot as plt
import copy
import py_model.hardio as hio
import graph_help as gh

NUM_TSTEPS = 18
NUM_CELLS = 2
COA = 0
CIL = 0
NUM_INT_STEPS = 10

py_out = read_save_file("py", NUM_TSTEPS, NUM_INT_STEPS, NUM_CELLS, CIL, COA)
rust_out = read_save_file("rust", NUM_TSTEPS, NUM_INT_STEPS, NUM_CELLS, CIL,
                          COA)

check_header_equality(py_out, rust_out)

py_dat = py_out["tsteps"]
rust_dat = rust_out["tsteps"]

py_data_dict_per_cell = gen_data_dict_per_cell_int_steps(py_out)
rust_data_dict_per_cell = gen_data_dict_per_cell_int_steps(rust_out)

py_data_dict_per_cell, rust_data_dict_per_cell = gh.prep_data(
    py_data_dict_per_cell, rust_data_dict_per_cell)

ALL_GROUP = gh.PlotDataGroup(py_data_dict_per_cell, rust_data_dict_per_cell,
                             gh.ALL_LABELS_NO_EXTRAS, gh.ALL_LABELS_NO_EXTRAS,
                             "all data")

POSITION_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        ["poly"] + gh.flatten(gh.EXTRA_POSITIONS),
        ["poly"] + gh.EXTRA_POSITIONS,
        "position related data",
    )

RAC_ACT_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        ["rgtp_forces", "rac_acts", "kgtps_rac", "x_coas", "coa_update"],
        ["rgtp_forces", "rac_acts", "kgtps_rac", "x_coas", "coa_update"],
        "Rac activity related data",
    )

RHO_ACT_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        ["rgtp_forces", "rho_acts", "kgtps_rho", "x_cils", "cil_update"],
        ["rgtp_forces", "rho_acts", "kgtps_rho", "x_cils", "cil_update"],
        "Rac activity related data",
    )

RAC_ACT_CONC_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        ["conc_rac_acts", "rac_act_net_fluxes"],
        ["conc_rac_acts", "rac_act_net_fluxes"],
        "Rac activity concentration related data",
    )

RAC_INACT_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        ["rac_inacts", "x_cils", "kdgtps_rac", "rho_acts"],
        [("rac_inacts", "rho_acts"), "x_cils", "kdgtps_rac"],
        "Rac inactivity related data",
    )

RHO_INACT_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        ["rho_inacts", "kdgtps_rho", "rac_acts", "x_coas"],
        [("rho_inacts", "rac_acts"), "kdgtps_rho", "x_coas"],
        "Rho inactivity related data",
    )

all_force_labels = ["sum_forces"] + hio.FORCES
FORCE_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        all_force_labels + gh.flatten(gh.EXTRA_FORCES),
        [tuple(all_force_labels), ("edge_forces", "edge_forces_minus")]
        + all_force_labels + gh.EXTRA_FORCES,
        "force related data",
    )

INTERACTIONS_GROUP = \
    gh.PlotDataGroup(
        py_data_dict_per_cell,
        rust_data_dict_per_cell,
        ["coa_update", "x_coas", "cil_update", "x_cils"],
        ["coa_update", "x_coas", "cil_update", "x_cils"],
        "COA and CIL related data",
    )

DATA_GROUPS = [POSITION_GROUP, FORCE_GROUP, RAC_ACT_GROUP, RHO_ACT_GROUP,
               RAC_INACT_GROUP, RHO_INACT_GROUP, INTERACTIONS_GROUP, ALL_GROUP]

PLOT_X_MIN = 0
PLOT_X_MAX = int(NUM_TSTEPS * NUM_INT_STEPS)
VERTEX_PLOT_TYPE = [n for n in range(16)] + ["all"]
CELL_PLOT_TYPE = [m for m in range(NUM_CELLS)] + ["all"]
VERT_PLOT_IX = 0
CELL_PLOT_IX = 0

DATA_GROUP_IX = 0
ACTIVE_DG = DATA_GROUPS[DATA_GROUP_IX]
NUM_LABEL_GROUPS = len(ACTIVE_DG.label_groups)
CURR_INNER_IX = 0


def paint(delta_vertex_plot, delta_data_plot, delta_max_plot,
          delta_cell_plot, delta_data_group_plot):
    global fig
    global ax
    global VERT_PLOT_IX
    global CURR_INNER_IX
    global PLOT_X_MAX
    global DATA_GROUP_IX
    global CELL_PLOT_IX
    global NUM_CELLS
    global DATA_GROUPS
    global ACTIVE_DG
    global NUM_LABEL_GROUPS
    ax.cla()

    if abs(delta_data_group_plot) > 0:
        DATA_GROUP_IX = (DATA_GROUP_IX + delta_data_group_plot) % len(
            DATA_GROUPS)
        ACTIVE_DG = DATA_GROUPS[DATA_GROUP_IX]
        NUM_LABEL_GROUPS = len(ACTIVE_DG.label_groups)
        CURR_INNER_IX = 0

    CURR_INNER_IX = (CURR_INNER_IX + delta_data_plot) % NUM_LABEL_GROUPS
    label_group = ACTIVE_DG.label_groups[CURR_INNER_IX]

    VERT_PLOT_IX = (VERT_PLOT_IX + delta_vertex_plot) % len(
        VERTEX_PLOT_TYPE)
    CELL_PLOT_IX = (CELL_PLOT_IX + delta_cell_plot) % len(CELL_PLOT_TYPE)

    vert = VERTEX_PLOT_TYPE[VERT_PLOT_IX]
    cell = CELL_PLOT_TYPE[CELL_PLOT_IX]

    if abs(delta_max_plot) > 0:
        PLOT_X_MAX = (PLOT_X_MAX + delta_max_plot) % (
                    NUM_TSTEPS * NUM_INT_STEPS)
        for g in DATA_GROUPS:
            g.recalc_ylims(PLOT_X_MAX)

    for m in range(NUM_CELLS):
        if m == cell or cell == "all":
            py_cell_data = ACTIVE_DG.py_dat[m]
            rust_cell_data = ACTIVE_DG.rust_dat[m]
            if type(label_group) != tuple:
                tupleized_label_group = (label_group,)
            else:
                tupleized_label_group = label_group

            for label in tupleized_label_group:
                color = gh.LABEL_COLOR_DICT[label]
                if len(py_cell_data[label].shape) == 1:
                    ax.plot(
                        rust_cell_data[label][:PLOT_X_MAX],
                        color=color, label=label)
                    ax.plot(
                        py_cell_data[label][:PLOT_X_MAX],
                        color=color,
                        linestyle="dashed", label=label)
                   #ax.set_ylim(ACTIVE_DG.grouped_ylims_dict[label_group])
                else:
                    for n in range(16):
                        if n == vert or vert == "all":
                            ax.plot(
                                rust_cell_data[label][:PLOT_X_MAX, n],
                                color=color, label=label)
                            ax.plot(
                                py_cell_data[label][:PLOT_X_MAX, n],
                                color=color,
                                linestyle="dashed", label=label)
                            #ax.set_ylim(ACTIVE_DG.grouped_ylims_dict[
                            # label_group])

    inter_tick_len = np.max([1, np.ceil(PLOT_X_MAX / 20)])
    xticks = np.arange(0, PLOT_X_MAX, inter_tick_len)[:20]

    ax.set_xticks(xticks)
    ax.grid(which="major", axis="x")
    ax.legend(loc="best")
    ax.set_title("{}\n{}, vert: {}, cell: {}".format(ACTIVE_DG.description,
                                                     label_group,
                                                     vert, cell))


def on_press(event):
    global VERT_PLOT_IX
    global CURR_INNER_IX
    global PLOT_X_MAX
    global DATA_GROUP_IX
    global CELL_PLOT_IX
    global fig
    print("pressed: {}".format(event.key))
    if event.key == "down":
        paint(-1, 0, 0, 0, 0)
    elif event.key == "up":
        paint(1, 0, 0, 0, 0)
    elif event.key == "left":
        paint(0, -1, 0, 0, 0)
    elif event.key == "right":
        paint(0, 1, 0, 0, 0)
    elif event.key == "z":
        paint(0, 0, -1, 0, 0)
    elif event.key == "x":
        paint(0, 0, 1, 0, 0)
    elif event.key == "c":
        paint(0, 0, -10, 0, 0)
    elif event.key == "v":
        paint(0, 0, 10, 0, 0)
    elif event.key == "b":
        paint(0, 0, -1000, 0, 0)
    elif event.key == "n":
        paint(0, 0, 1000, 0, 0)
    elif event.key == "pagedown":
        paint(0, 0, 0, 0, -1)
    elif event.key == "pageup":
        paint(0, 0, 0, 0, 1)
    elif event.key == "r":
        VERT_PLOT_IX = 0
        CURR_INNER_IX = 0
        #PLOT_X_MIN = 0
        PLOT_X_MAX = NUM_TSTEPS * NUM_INT_STEPS
        CELL_PLOT_IX = 0
        DATA_GROUP_IX = 0
        paint(0, 0, 0, 0, 0)


fig, ax = plt.subplots()
paint(0, 0, 0, 0, 0)
fig.canvas.mpl_connect('key_press_event', on_press)
