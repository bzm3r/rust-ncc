use crate::NVERTS;

/// `CcVvDat` allows storage of (cell, vertex)-(cell, vertex) data
/// indexed by `(ci, vi, oci, ovi)`, (`ci` for "cell index", `oci` for
/// "other cell index", `vi` for "vertex index", and `ovi` for "other
/// vertex index"). Data is stored in chunks relative to one cell, and
/// another (guaranteed not to be itself). All data stored must have
/// the same type. However, the structure is generic over different
/// data types as long as they implement `Copy` and `Default`.
#[allow(unused)]
pub struct CcVvDat<T: Copy> {
    /// Number of cells this structure is designed to store data for.
    num_cells: usize,
    vv_stride: usize,
    cc_stride: usize,
    /// Vector in which data is stored linearly.
    dat: Vec<T>,
    /// Value to return when no data is set, or when an invalid
    /// index of the form `(ci, vi, ci, ovi)`.
    undefined: T,
}

#[allow(unused)]
impl<T: Copy> CcVvDat<T> {
    /// Generate an empty `CcVvDat` structure.
    pub fn empty(num_cells: usize, undefined: T) -> CcVvDat<T> {
        // number of possible vertex-vertex pairs.
        let vv_stride = NVERTS * NVERTS;
        // number of cells a cell can interact with (ignores self)
        let cc_stride = (num_cells - 1) * vv_stride;
        CcVvDat {
            num_cells,
            vv_stride,
            cc_stride,
            dat: vec![undefined; num_cells * cc_stride],
            undefined,
        }
    }

    /// Calculate the index of the element in `dat` which corresponds
    /// `(ci, vi, oci, ovi)`.
    pub fn calc_ix(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        ovi: usize,
    ) -> usize {
        let oci = if ci < oci { oci - 1 } else { oci };
        ci * self.cc_stride + oci * self.vv_stride + vi * NVERTS + ovi
    }

    /// Set data into the element indexed by `(ci, vi, oci, ovi)`.
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

    /// Set data into the element indexed by `(ci, vi, oci, ovi)`.
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

    /// Get data from the element indexed by `(ci, vi, oci, ovi)`.
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

/// `CvCvDat` allows storage of (cell, vertex)-(cell, vertex) data
/// indexed by `(ci, vi, oci, ovi)`, (`ci` for "cell index", `oci` for
/// "other cell index", `vi` for "vertex index", and `ovi` for "other
/// vertex index"). Data is stored in chunks relative to one cell, and
/// a particular vertex on it. All data stored must have
/// the same type. However, the structure is generic over different
/// data types as long as they implement `Copy` and `Default`.
#[derive(Clone)]
pub struct CvCvDat<T: Copy> {
    pub num_cells: usize,
    c_stride: usize,
    cv_stride: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy> CvCvDat<T> {
    /// Generate an empty `CvCvDat` structure.
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

    pub fn calc_ix(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        ovi: usize,
    ) -> usize {
        let oci = if ci < oci { oci - 1 } else { oci };
        ci * self.c_stride + vi * self.cv_stride + oci * NVERTS + ovi
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
