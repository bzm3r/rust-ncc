import numpy as np


def p2ds_to_numpy(p2ds):
    vs = []
    for p2d in p2ds:
        vs.append([p2d['x'], p2d['y']])
    return np.array(vs)


def extract_p2ds_from_data(address, data):
    dat_per_cell_per_tstep = []
    for cell_recs in data:
        dat_per_cell = []
        for cell_rec in cell_recs:
            for addr in address:
                cell_rec = cell_rec[addr]
            dat_per_cell.append(p2ds_to_numpy(cell_rec))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_scalars_from_data(address, data):
    dat_per_cell_per_tstep = []
    for cell_recs in data:
        dat_per_cell = []
        for cell_rec in cell_recs:
            for addr in address:
                cell_rec = cell_rec[addr]
            dat_per_cell.append(cell_rec)
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_p2ds_from_cell_states(state_key, dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['states']:
            dat_per_cell.append(p2ds_to_numpy(cell_rec[state_key][dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_p2ds_from_interactions(dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['interactions']:
            dat_per_cell.append(p2ds_to_numpy(cell_rec[dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_scalars(state_key, dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['states']:
            dat_per_cell.append(np.array(cell_rec[state_key][dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)
