use crate::interactions::utils::{sort_ixs4d, sym_sum};
use crate::NVERTS;

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
        let sym_sum_n = sym_sum(num_cells);
        let vv_stride = NVERTS * NVERTS;
        SymCcVvDat {
            num_cells,
            sym_sum_n,
            vv_stride,
            dat: vec![undefined; sym_sum_n * vv_stride],
            undefined,
        }
    }

    pub fn calc_ix(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        ovi: usize,
    ) -> usize {
        let (ci, vi, oci, ovi) = sort_ixs4d(ci, vi, oci, ovi);
        let k = self.sym_sum_n - sym_sum(self.num_cells - ci);
        let ixc = k + (oci - ci - 1);
        ixc * self.vv_stride + vi * NVERTS + ovi
    }

    pub fn set(
        &mut self,
        ci: usize,
        vi: usize,
        oci: usize,
        ovi: usize,
        x: T,
    ) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = x;
        }
    }

    #[allow(unused)]
    pub fn reset(
        &mut self,
        ci: usize,
        vi: usize,
        oci: usize,
        ovi: usize,
    ) {
        if ci != oci {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix] = self.undefined;
        }
    }

    pub fn get(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        ovi: usize,
    ) -> T {
        if ci == oci {
            self.undefined
        } else {
            let ix = self.calc_ix(ci, vi, oci, ovi);
            self.dat[ix]
        }
    }
}
