use crate::cell::ModelCell;
use crate::math::geometry::Bbox;
use crate::math::p2d::P2D;
use crate::math::{min_f32s, min_f32s_ix};
use crate::parameters::Parameters;
use crate::utils::circ_ix_plus;
use crate::world::Interactions;
use crate::NVERTS;

pub struct Mat<T: Copy + Default> {
    n_r: usize,
    n_c: usize,
    dat: Vec<T>,
}

impl<T: Copy + Default> Mat<T> {
    pub fn new(n_r: usize, n_c: usize, default: T) -> Mat<T> {
        Mat {
            n_r,
            n_c,
            dat: vec![default; n_r * n_c],
        }
    }

    fn calc_ix(&self, i: usize, j: usize) -> usize {
        i * self.n_r + j
    }

    pub fn store(&mut self, i: usize, j: usize, x: T) {
        self.dat[self.calc_ix(i, j)] = x;
    }

    pub fn get(&self, i: usize, j: usize) -> T {
        self.dat[self.calc_ix(i, j)]
    }
}

pub struct SymMat<T: Copy + Default> {
    n: usize,
    dat: Vec<T>,
}

impl<T: Copy + Default> SymMat<T> {
    pub fn new(n: usize, default: T) -> SymMat<T> {
        SymMat {
            n,
            dat: vec![default; n + (n * (n - 1) / 2)],
        }
    }

    fn calc_ix(&self, i: usize, j: usize) -> usize {
        let (i, j) = if i < j { (i, j) } else { (j, i) };
        i * self.n + j
    }

    pub fn set(&mut self, i: usize, j: usize, x: T) {
        self.dat[self.calc_ix(i, j)] = x;
    }

    pub fn get(&self, i: usize, j: usize) -> T {
        self.dat[self.calc_ix(i, j)]
    }
}

/// Stores distance between vertices of cells.
pub struct DistMat(SymMat<f32>);

impl DistMat {
    /// Get the distance between vertex `vi` of cell `ci` and vertex `ovi` of cell `oci`.
    pub fn get(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> f32 {
        self.0.get(ci * NVERTS + vi, oci * NVERTS + ovi)
    }

    /// Get `(ix, ix2, d)`, where `ix` and `ix2` are the vertices on `oci` defining the edge
    /// containing the point closest to the vertex `vi` of `ci`.
    pub fn get_close_point_info(&self, ci: usize, vi: usize, oci: usize) -> (usize, usize, f32) {
        let dists = (0..NVERTS)
            .map(|ovi| self.get(ci, vi, oci, ovi))
            .collect::<Vec<f32>>();
        let (ix, d) = min_f32s_ix(&dists);
        (ix, circ_ix_plus(ix, NVERTS), d)
    }
}

/// Calculate a symmetric matrix whose `(i, j)` element is a boolean representing whether cells `i`
/// and `j` are in contact range.
pub fn calc_contact_mat(cells: &[ModelCell], contact_range: f32) -> SymMat<bool> {
    let num_cells: usize = cells.len();
    let bboxes = cells
        .iter()
        .map(|c| Bbox::calc(&c.state.vertex_coords).expand_by(contact_range))
        .collect::<Vec<Bbox>>();
    let mut contact_mat = SymMat::new(num_cells, false);
    for ci in 0..num_cells {
        let this_bbox = &bboxes[ci];
        for oci in (ci + 1)..num_cells {
            contact_mat.set(ci, oci, bboxes[oci].intersects(this_bbox));
        }
    }
    contact_mat
}

/// Calculate distance between vertices of cells. Distance is `f32::INFINITY` the cells are not in
/// contact, as defined by the supplied contact matrix.
pub fn calc_dist_mat(
    cells: &[ModelCell],
    contact_mat: &SymMat<bool>,
    contact_range: f32,
) -> DistMat {
    let num_cells: usize = cells.len();
    let mut dist_mat = SymMat::new(num_cells * NVERTS, (0, f32::INFINITY));
    for ci in 0..num_cells {
        let this_vcs = &cells[ci].state.vertex_coords;
        for oci in (ci + 1)..num_cells {
            if contact_mat.get(ci, oci) {
                let other_vcs = &cells[oci].state.vertex_coords;
                for vi in 0..(NVERTS) {
                    for ovi in 0..(NVERTS) {
                        let ovi2 = circ_ix_plus(ovi, NVERTS);
                        let d = calc_dist_of_point_to_seg(
                            &this_vcs[vi],
                            &other_vcs[ovi],
                            &other_vcs[ovi2],
                        );
                    }
                }
            }
        }
    }
    DistMat(dist_mat)
}

#[derive(Clone, Default)]
pub struct CloseCellInfo {
    num_close_cells: usize,
    closest_cells: [usize; 9],
}

impl CloseCellInfo {
    pub fn insert(&mut self, ix: usize) {
        if self.num_close_cells < 9 {
            self.closest_cells[self.num_close_cells] = ix;
            self.num_close_cells += 1;
        }
    }
}

impl Interactions {
    pub fn generate(cells: &[ModelCell], parameters: &Parameters) -> Vec<Interactions> {
        let contact_mat = calc_contact_mat(cells, parameters.close_criterion);
        let dist_mat = calc_dist_mat(cells, &contact_mat);
        let num_cells = cells.len();
        let mut closest_cells = vec![CloseCellInfo::default(); num_cells];
        let mut interactions = vec![Interactions::default(); num_cells];
        for ci in 0..num_cells {
            let this_vcs = &cells[ci].state.vertex_coords;
            for oci in 0..num_cells {
                if contact_mat.get(ci, oci) {
                    closest_cells[ci].insert(oci);
                    let other_vcs = &cells[ci].state.vertex_coords;
                    for vi in 0..NVERTS {
                        let (ovi, ovi2, dist) = dist_mat.get_closest_vertex(ci, vi, oci);
                    }
                }
            }
        }
        interactions
    }
}
