# -*- coding: utf-8 -*-
"""
Created on Sun Jun 07 17:47:39 2015

@author: Brian
"""

import numpy as np


# ==============================================================
def calculate_normalized_randomization_factors(size):
    rfs = np.random.random(size)
    return rfs / np.sum(rfs)


# ===============================================================


def is_numeric(s):
    try:
        float(s)
        return True

    except ValueError:
        return False

# @nb.jit(nopython=True)


def make_verts_array_given_xs_and_ys(xs, ys):
    verts = np.empty((16, 2), dtype=np.float64)

    for i in range(16):
        vert = verts[i]
        vert[0] = xs[i]
        vert[1] = ys[i]

    return verts
