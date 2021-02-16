# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 16:00:16 2016

@author: Brian Merchant
"""

import parameters
import environment
import numpy as np


def rust_comparison_test(
        raw_params,
        out_dir,
        name,
        log_level,
        box_width=2,
        box_height=1,
):
    params = parameters.refine_raw_params(raw_params)
    cell_d = 2 * params["cell_r"]

    box_height = box_height * cell_d
    box_width = box_width * cell_d
    box_x_offset = 0.0
    box_y_offset = 0.0

    cell_group_bbox = np.array(
        [
            box_x_offset,
            box_x_offset + box_width,
            box_y_offset,
            box_height + box_y_offset,
            ]
    )
    params["cell_group_bbox"] = cell_group_bbox

    an_environment = environment.Environment(out_dir, name,
                                             log_level, params)

    an_environment.execute_system_dynamics()

    print("Done.")

    return params
