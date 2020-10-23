use crate::NVERTS;

fn sym_sum_to(n: usize) -> usize {
    n * (n - 1) / 2
}

fn fix_oci(ci: usize, oci: usize) -> usize {
    if ci < oci {
        oci - 1
    } else {
        oci
    }
}

fn sort_ixs(ci: usize, oci: usize) -> (usize, usize) {
    if oci < ci {
        (oci, ci)
    } else {
        (ci, oci)
    }
}

fn sort_ixs4d(ci: usize, vi: usize, oci: usize, ovi: usize) -> (usize, usize, usize, usize) {
    if oci < ci {
        (oci, ovi, ci, vi)
    } else {
        (ci, vi, oci, ovi)
    }
}

#[cfg(debug_assertions)]
pub fn check_indices4d(num_cells: usize, ci: usize, vi: usize, oci: usize, ovi: usize) {
    let oci = fix_oci(ci, oci);

    if ci > num_cells {
        panic!("{} cells tracked, received ci: {}", num_cells, ci);
    }

    if vi > NVERTS as usize {
        panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
    }

    if oci > num_cells - 1 {
        panic!(
            "{} cells tracked, received oci: {} (max: {})",
            num_cells,
            oci,
            num_cells - 1
        );
    }

    if ovi > NVERTS as usize {
        panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
    }
}

#[cfg(debug_assertions)]
pub fn check_sym_indices4d(num_cells: usize, ci: usize, vi: usize, oci: usize, ovi: usize) {
    let (ci, vi, oci, ovi) = sort_ixs4d(ci, vi, oci, ovi);

    if ci > num_cells {
        panic!("{} cells tracked, received ci: {}", num_cells, ci);
    }

    if vi > NVERTS as usize {
        panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
    }

    if oci > num_cells {
        panic!(
            "{} cells tracked, received symmetric oci: {}",
            num_cells, oci
        );
    }

    if ovi > NVERTS as usize {
        panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
    }
}

#[derive(Clone, Default)]
pub struct SymCcDat<T: Copy + Default> {
    pub num_cells: usize,
    sym_sum_nc: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy + Default> SymCcDat<T> {
    pub fn new(num_cells: usize, undefined: T) -> SymCcDat<T> {
        let sym_sum_nc = sym_sum_to(num_cells);
        SymCcDat {
            num_cells,
            sym_sum_nc,
            dat: vec![undefined; sym_sum_nc],
            undefined,
        }
    }

    fn calc_ix(&self, ci: usize, oci: usize) -> usize {
        let (ci, oci) = sort_ixs(ci, oci);
        let k = sym_sum_to(self.num_cells) - sym_sum_to(self.num_cells - ci);
        k + (ci - oci - 1)
    }

    pub fn set(&mut self, ci: usize, oci: usize, x: T) {
        if ci != oci {
            let ix = self.calc_ix(ci, oci);
            self.dat[ix] = x;
        }
    }

    pub fn get(&self, ci: usize, oci: usize) -> T {
        if ci == oci {
            self.undefined
        } else {
            let ix = self.calc_ix(ci, oci);
            self.dat[ix]
        }
    }
}

#[derive(Clone, Default)]
pub struct CcDat<T: Copy + Default> {
    pub num_cells: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy + Default> CcDat<T> {
    pub fn new(num_cells: usize, undefined: T) -> CcDat<T> {
        CcDat {
            num_cells,
            dat: vec![undefined; num_cells * (num_cells - 1)],
            undefined,
        }
    }

    fn calc_ix(&self, ci: usize, oci: usize) -> usize {
        let oci = fix_oci(ci, oci);
        ci * (self.num_cells - 1) + oci
    }

    pub fn set(&mut self, ci: usize, oci: usize, x: T) {
        if ci != oci {
            let ix = self.calc_ix(ci, oci);
            self.dat[ix] = x;
        }
    }

    pub fn get(&self, ci: usize, oci: usize) -> T {
        if ci == oci {
            self.undefined
        } else {
            let ix = self.calc_ix(ci, oci);
            self.dat[ix]
        }
    }
}

/// Matrix to store inter-vertex data in `[cell][cell][vertex][vertex]` format.
pub struct CcVvDat<T: Copy> {
    num_cells: usize,
    vv_stride: usize,
    cc_stride: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy> CcVvDat<T> {
    pub fn empty(num_cells: usize, undefined: T) -> CcVvDat<T> {
        let vv_stride = NVERTS * NVERTS;
        let cc_stride = (num_cells - 1) * vv_stride;
        CcVvDat {
            num_cells,
            vv_stride,
            cc_stride,
            dat: vec![undefined; num_cells * cc_stride],
            undefined,
        }
    }

    pub fn calc_ix(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> usize {
        let oci = fix_oci(ci, oci);
        ci * self.cc_stride + oci * self.vv_stride + vi * NVERTS + ovi
    }

    pub fn set(&mut self, ci: usize, vi: usize, oci: usize, ovi: usize, x: T) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = x;
        }
    }

    pub fn reset(&mut self, ci: usize, vi: usize, oci: usize, ovi: usize) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = self.undefined;
        }
    }

    pub fn get(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> T {
        if ci == oci {
            self.undefined
        } else {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix]
        }
    }
}

/// Matrix to store inter-vertex data.
#[derive(Clone)]
pub struct SymCcVvDat<T: Copy> {
    num_cells: usize,
    sym_sum_n: usize,
    vv_stride: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy> SymCcVvDat<T> {
    pub fn empty(num_cells: usize, undefined: T) -> SymCcVvDat<T> {
        let sym_sum_n = sym_sum_to(num_cells);
        let vv_stride = NVERTS * NVERTS;
        SymCcVvDat {
            num_cells,
            sym_sum_n,
            vv_stride,
            dat: vec![undefined; sym_sum_n * vv_stride],
            undefined,
        }
    }

    pub fn calc_ix(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> usize {
        let (ci, vi, oci, ovi) = sort_ixs4d(ci, vi, oci, ovi);
        let k = self.sym_sum_n - sym_sum_to(self.num_cells - ci);
        let ixc = k + (oci - ci - 1);
        ixc * self.vv_stride + vi * NVERTS + ovi
    }

    pub fn set(&mut self, ci: usize, vi: usize, oci: usize, ovi: usize, x: T) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = x;
        }
    }

    pub fn reset(&mut self, ci: usize, vi: usize, oci: usize, ovi: usize) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = self.undefined;
        }
    }

    pub fn get(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> T {
        if ci == oci {
            self.undefined
        } else {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix]
        }
    }
}

/// Matrix to store inter-vertex data in `[cell][cell][vertex][vertex]` format.
#[derive(Clone)]
pub struct CvCvDat<T: Copy> {
    pub num_cells: usize,
    c_stride: usize,
    cv_stride: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy> CvCvDat<T> {
    pub fn empty(num_cells: usize, undefined: T) -> CvCvDat<T> {
        let cv_stride = (num_cells - 1) * NVERTS;
        let c_stride = NVERTS * cv_stride;
        CvCvDat {
            num_cells,
            c_stride,
            cv_stride,
            dat: vec![undefined; num_cells * c_stride],
            undefined,
        }
    }

    pub fn calc_ix(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> usize {
        let oci = fix_oci(ci, oci);
        ci * self.c_stride + vi * self.cv_stride + oci * NVERTS + ovi
    }

    pub fn set(&mut self, ci: usize, vi: usize, oci: usize, ovi: usize, x: T) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = x;
        }
    }

    pub fn reset(&mut self, ci: usize, vi: usize, oci: usize, ovi: usize) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = self.undefined;
        }
    }

    pub fn get(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> T {
        if ci == oci {
            self.undefined
        } else {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix]
        }
    }
}
