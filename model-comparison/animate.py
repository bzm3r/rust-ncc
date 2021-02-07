from retrieve import *
import numpy as np
import matplotlib.pyplot as plt

# ===========================================================================

NUM_TSTEPS = 1800
NUM_CELLS = 2
COA = 24
CIL = 60
NUM_INT_STEPS = 10

py_out = read_save_file("py", NUM_TSTEPS, NUM_INT_STEPS, NUM_CELLS, CIL, COA)
rust_out = read_save_file("rust", NUM_TSTEPS, NUM_INT_STEPS, NUM_CELLS, CIL,
                          COA)
check_header_equality(py_out, rust_out)

py_data_dict_per_cell = gen_data_dict_per_cell_tsteps(py_out)
rust_data_dict_per_cell = gen_data_dict_per_cell_tsteps(rust_out)


# ===========================================================================

def gen_animation_material(data_dict_per_cell):
    poly_per_cell_per_tstep = np.array([[data_dict_per_cell[c]["poly"][t]
                                         for c in range(NUM_CELLS)]
                                        for t in range(NUM_TSTEPS)])
    uivs_per_cell_per_tstep = np.array([[data_dict_per_cell[c]["uivs"][t]
                                         for c in range(NUM_CELLS)]
                                        for t in range(NUM_TSTEPS)])
    uovs_per_cell_per_step = -1 * uivs_per_cell_per_tstep
    rac_acts_per_cell_per_tstep = np.array([[data_dict_per_cell[c][
                                                 "rac_acts"][t]
                                             for c in range(NUM_CELLS)]
                                            for t in range(NUM_TSTEPS)])
    rac_act_arrows_per_cell_per_tstep = \
        50 * rac_acts_per_cell_per_tstep[:, :, :, np.newaxis] * \
        uovs_per_cell_per_step
    rho_acts_per_cell_per_tstep = np.array([[data_dict_per_cell[c][
                                                 "rho_acts"][t]
                                             for c in range(NUM_CELLS)]
                                            for t in range(NUM_TSTEPS)])
    rho_act_arrows_per_cell_per_tstep = \
        50 * rho_acts_per_cell_per_tstep[:, :, :, np.newaxis] * \
        uivs_per_cell_per_tstep

    return poly_per_cell_per_tstep, uivs_per_cell_per_tstep, \
           uovs_per_cell_per_step, rac_acts_per_cell_per_tstep, \
           rac_act_arrows_per_cell_per_tstep, rho_act_arrows_per_cell_per_tstep


# ===========================================================================

py_ani_mats = gen_animation_material(py_data_dict_per_cell)
rust_ani_mats = gen_animation_material(rust_data_dict_per_cell)


# ===========================================================================

def paint_cell(ax, ani_mats, poly_ls):
    poly_per_cell_per_tstep, uivs_per_cell_per_tstep, \
    uovs_per_cell_per_step, rac_acts_per_cell_per_tstep, \
    rac_act_arrows_per_cell_per_tstep, rho_act_arrows_per_cell_per_tstep = \
        ani_mats
    for (ci, poly) in enumerate(poly_per_cell_per_tstep[tstep]):
        if ci == 0:
            poly_color = "k"
        else:
            poly_color = "g"

        for vix in range(16):
            ax.plot([poly[vix, 0], poly[(vix + 1) % 16, 0]],
                    [poly[vix, 1], poly[(vix + 1) % 16, 1]],
                    color=poly_color, linestyle=poly_ls)
            ax.annotate(str(vix), (poly[vix, 0], poly[vix, 1]))

    for poly, rac_act_arrows in zip(
            poly_per_cell_per_tstep[tstep],
            rac_act_arrows_per_cell_per_tstep[tstep]
    ):
        for p, rac_arrow in zip(poly, rac_act_arrows):
            ax.arrow(p[0], p[1], 1 * rac_arrow[0], 1 * rac_arrow[1],
                     color="b",
                     length_includes_head=True, head_width=0.0)

    for poly, rho_act_arrows in zip(poly_per_cell_per_tstep[tstep],
                                    rho_act_arrows_per_cell_per_tstep[
                                        tstep]):
        for p, rho_arrow in zip(poly, rho_act_arrows):
            ax.arrow(p[0], p[1], 1 * rho_arrow[0], 1 * rho_arrow[1],
                     color="r",
                     length_includes_head=True, head_width=0.0)


# ===========================================================================

def paint(delta):
    global fig
    global ax
    global tstep
    global NUM_TSTEPS
    global py_ani_mats
    global rust_ani_mats
    ax.cla()
    ax.set_aspect('equal')
    ax.set_xlim([-40, 200])
    ax.set_ylim([-40, 200])

    paint_cell(ax, rust_ani_mats, "-")
    paint_cell(ax, py_ani_mats, ":")

    ax.set_title("frame {}".format(tstep))
    tstep = (tstep + delta) % NUM_TSTEPS
    plt.show()


# ===========================================================================

def on_press(event):
    global fig
    if event.key == 'x':
        paint(1)
    elif event.key == 'z':
        paint(-1)
    if event.key == 'c':
        paint(-5)
    elif event.key == 'v':
        paint(5)
    elif event.key == 'n':
        paint(-10)
    elif event.key == 'm':
        paint(10)
    fig.canvas.draw()


# ===========================================================================

tstep = 0
fig, ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', on_press)
paint(0)
