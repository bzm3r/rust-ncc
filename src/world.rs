// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::cell::{Cell, CellState};
use std::fs::File;
use std::ops::Range;
use std::path::Path;
use crate::experiment::{Experiment, GroupLayout};
use crate::parameters::{UserParams, CharQuants};
use crate::math::P2D;
use crate::quantity::Length;

struct WorldState {
    curr_tstep: usize,
    cells: Vec<Cell>,
}

impl WorldState {
    fn simulate(&mut self, final_tstep: usize) {
        while self.curr_tstep < final_tstep {
            #[cfg(debug_assertions)]
                println!("tstep: {}\n==========", self.curr_tstep);
            for cell in &mut self.cells {
                #[cfg(debug_assertions)]
                    println!("executing cell: {}", cell.ix);
                cell.simulate();
                #[cfg(debug_assertions)]
                    println!("executing cell: {:?}", cell.state);
            }
            self.curr_tstep += 1;
            #[cfg(debug_assertions)]
                println!("--------");
        }
    }
}

// struct WorldHistory {
//     fd: File,
//     live: Vec<WorldHistory>,
//     max_live: usize,
// }

// impl WorldHistory {
//     fn new(_path: &Path) -> WorldHistory {
//         unimplemented!()
//     }
//
//     fn send(&mut self, _state: WorldState) {
//         unimplemented!()
//     }
//
//     fn retrieve_state(&self, _tstep_range: Range<usize>) -> Vec<WorldState> {
//         unimplemented!()
//     }
//
//     fn retrieve_cell_state(&self, _tstep_range: Range<usize>, _cix: usize) -> Vec<CellState> {
//         unimplemented!()
//     }
// }

// impl Drop for WorldHistory {
//     fn drop(&mut self) {
//         unimplemented!()
//     }
// }

pub(crate) struct World {
    //history: WorldHistory,
    tstep_in_secs: f32,
    state: WorldState,
}


fn gen_cell_centroids(layout: &GroupLayout, num_cells: u32, cell_diam: Length, char_quants: &CharQuants) -> Result<Vec<P2D>, String> {
    let cell_diam = char_quants.normalize(&cell_diam);
    let group_centroid = P2D {
        x: char_quants.normalize(&layout.centroid[0]),
        y: char_quants.normalize(&layout.centroid[1]),
    };
    if layout.width * layout.height >= num_cells {
        let nr = layout.height / num_cells;
        let nc = layout.width / num_cells;
        let mut r = vec![];
        let first_cell_centroid = {
            let delta = P2D {
                x: ((layout.width - 1) as f32) * 0.5 * cell_diam,
                y: ((layout.height - 1) as f32) * 0.5 * cell_diam,
            };
            group_centroid - delta
        };
        let cd = P2D {
            x: 0.5 * cell_diam,
            y: 0.5 * cell_diam,
        };
        for ix in 0..num_cells {
            let row = (ix / nr) as f32;
            let col = (ix as f32) - row;
            r.push(first_cell_centroid + cd.scalar_mulx(row) + cd.scalar_muly(col));
        }
        Ok(r)
    } else {
        Err(format!("Cell group layout area ({}x{}) not large enough to fit {} cells.", layout.width, layout.height, num_cells))
    }
}

impl World {
    pub(crate) fn new(experiment: Experiment) -> World {
        let mut num_cells = 0;
        let mut cells = vec![];
        let char_quants = CharQuants::default();
        for (usize, cg) in experiment.cell_groups.iter().enumerate() {
            let user_params = UserParams::default().apply_overrides(&cg.parameter_overrides);
            let params = user_params.gen_params(&char_quants);
            let cell_centroids = gen_cell_centroids(&cg.layout, cg.num_cells, user_params.cell_diam, &char_quants).unwrap();
            for cc in cell_centroids.into_iter() {
                cells.push(Cell::new(num_cells, usize, params.clone(), cc));
                num_cells += 1;
            }
        }
        World {
            tstep_in_secs: char_quants.time(),
            state: WorldState {
                curr_tstep: 0,
                cells
            }
        }
    }

    // fn save_history(&self) {
    //     unimplemented!()
    // }

    pub(crate) fn simulate(&mut self, final_tpoint: f32) {
        let final_tstep = (final_tpoint / self.tstep_in_secs).ceil() as usize;
        self.state.simulate(final_tstep);
    }
}
