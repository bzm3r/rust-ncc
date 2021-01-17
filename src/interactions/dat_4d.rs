use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

// /// `CcVvDat` allows storage of (cell, vertex)-(cell, vertex) data
// /// indexed by `(ci, vi, oci, ovi)`, (`ci` for "cell index", `oci` for
// /// "other cell index", `vi` for "vertex index", and `ovi` for "other
// /// vertex index"). Data is stored in chunks relative to one cell, and
// /// another (guaranteed not to be itself). All data stored must have
// /// the same type. However, the structure is generic over different
// /// data types as long as they implement `Copy` and `Default`.
// pub struct CcVvDat<T: Copy> {
//     /// Number of cells this structure is designed to store data for.
//     num_cells: u32,
//     vv_stride: u32,
//     cc_stride: u32,
//     /// Vector in which data is stored linearly.
//     dat: Vec<T>,
//     /// Value to return when no data is set, or when an invalid
//     /// index of the form `(ci, vi, ci, ovi)`.
//     undefined: T,
// }
//
// impl<T: Copy> CcVvDat<T> {
//     /// Generate an empty `CcVvDat` structure.
//     pub fn empty(num_cells: u32, undefined: T) -> CcVvDat<T> {
//         // number of possible vertex-vertex pairs.
//         let vv_stride = NVERTS * NVERTS;
//         // number of cells a cell can interact with (ignores self)
//         let cc_stride = (num_cells - 1) * vv_stride;
//         CcVvDat {
//             num_cells,
//             vv_stride,
//             cc_stride,
//             dat: vec![undefined; num_cells * cc_stride],
//             undefined,
//         }
//     }
//
//     /// Calculate the index of the element in `dat` which corresponds
//     /// `(ci, vi, oci, ovi)`.
//     pub fn calc_ix(
//         &self,
//         ci: u32,
//         vi: u32,
//         oci: u32,
//         ovi: u32,
//     ) -> u32 {
//         let oci = if ci < oci { oci - 1 } else { oci };
//         ci * self.cc_stride + oci * self.vv_stride + vi * NVERTS + ovi
//     }
//
//     /// Set data into the element indexed by `(ci, vi, oci, ovi)`.
//     pub fn set(
//         &mut self,
//         ci: u32,
//         vi: u32,
//         oci: u32,
//         ovi: u32,
//         x: T,
//     ) {
//         if ci != oci {
//             let ix = self.calc_ix(ci, vi, oci, ovi);
//             self.dat[ix] = x;
//         }
//     }
//
//     /// Set data into the element indexed by `(ci, vi, oci, ovi)`.
//     pub fn reset(
//         &mut self,
//         ci: u32,
//         vi: u32,
//         oci: u32,
//         ovi: u32,
//     ) {
//         if ci != oci {
//             let ix = self.calc_ix(ci, vi, oci, ovi);
//             self.dat[ix] = self.undefined;
//         }
//     }
//
//     /// Get data from the element indexed by `(ci, vi, oci, ovi)`.
//     pub fn get(
//         &self,
//         ci: u32,
//         vi: u32,
//         oci: u32,
//         ovi: u32,
//     ) -> T {
//         if ci == oci {
//             self.undefined
//         } else {
//             let ix = self.calc_ix(ci, vi, oci, ovi);
//             self.dat[ix]
//         }
//     }
// }

/// `CvCvDat` allows storage of (cell, vertex)-(cell, vertex) data
/// indexed by `(ci, vi, oci, ovi)`, (`ci` for "cell index", `oci` for
/// "other cell index", `vi` for "vertex index", and `ovi` for "other
/// vertex index"). Data is stored in chunks relative to one cell, and
/// a particular vertex on it. All data stored must have
/// the same type. However, the structure is generic over different
/// data types as long as they implement `Copy` and `Default`.
#[derive(Clone, Serialize, Deserialize, PartialEq, Default, Debug)]
pub struct CvCvDat<T: Copy + Default + Debug + PartialEq> {
    pub num_cells: usize,
    c_stride: usize,
    cv_stride: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy + Debug + Default + PartialEq> CvCvDat<T> {
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

    // pub fn calc_ix_range(
    //     &self,
    //     ci: u32,
    //     vi: u32,
    //     oci: u32,
    // ) -> (u32, u32) {
    //     let oci = if ci < oci { oci - 1 } else { oci };
    //     let begin =
    //         ci * self.c_stride + vi * self.cv_stride + oci * NVERTS;
    //     (begin, begin + NVERTS)
    // }
    //
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

    // pub fn reset(
    //     &mut self,
    //     ci: u32,
    //     vi: u32,
    //     oci: u32,
    //     ovi: u32,
    // ) {
    //     if ci != oci {
    //         let ix = self.calc_ix(ci, vi, oci, ovi);
    //         self.dat[ix] = self.undefined;
    //     }
    // }

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

    // pub fn get_per_other_vertex(
    //     &self,
    //     ci: u32,
    //     vi: u32,
    //     oci: u32,
    // ) -> [T; NVERTS] {
    //     let mut r = [self.undefined; NVERTS];
    //     if ci != oci {
    //         let (begin, end) = self.calc_ix_range(ci, vi, oci);
    //         r.copy_from_slice(&self.dat[begin..end]);
    //     }
    //     r
    // }
}
