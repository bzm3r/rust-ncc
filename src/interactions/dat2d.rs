/// `CcDat` allows storage of cell-cell data indexed by `(ci, oci)`,
/// (`ci` for "cell index") and (`oci` for "other cell index"). All
/// data stored must have the same type. However, the structure is
/// generic over different data types as long as they implement `Copy`
/// and `Default`.
#[derive(Clone, Default)]
pub struct CcDat<T: Copy + Default> {
    pub num_cells: usize,
    dat: Vec<T>,
    undefined: T,
}

#[allow(unused)]
impl<T: Copy + Default> CcDat<T> {
    /// Create a new `CcDat` structure that can hold data for
    /// `num_cells` cells. `undefined` is the value stored when no
    /// data has been set for a particular element. Note that the
    /// structure holds `num_cells * (num_cells - 1)` elements. This
    /// is because we assume that we do not need to store data for
    /// `(ci, ci)` We assume instead that we always have `(ci, oci)`
    /// such that `ci != oci`.
    pub fn new(num_cells: usize, undefined: T) -> CcDat<T> {
        CcDat {
            num_cells,
            dat: vec![undefined; num_cells * (num_cells - 1)],
            undefined,
        }
    }

    /// Calculate the index at which data for `(ci, oci)` is stored.
    fn calc_ix(&self, ci: usize, oci: usize) -> usize {
        // shift `oci` to account to ignore "diagonal" elements
        if ci < oci {
            oci - 1
        } else {
            oci
        };
        // multiply by `self.num_cells - 1` to ignore diagonals
        ci * (self.num_cells - 1) + oci
    }

    /// Set `x` as data for index `(ci, oci)`.
    /// Setting does nothing if `ci == oci`.
    pub fn set(&mut self, ci: usize, oci: usize, x: T) {
        if ci != oci {
            let ix = self.calc_ix(ci, oci);
            self.dat[ix] = x;
        }
    }

    /// Get the data stored at index `(ci, oci)`. If `ci == oci`,
    /// this will return the undefined value specified during
    /// creation of `CcDat`.
    pub fn get(&self, ci: usize, oci: usize) -> T {
        if ci == oci {
            self.undefined
        } else {
            let ix = self.calc_ix(ci, oci);
            self.dat[ix]
        }
    }
}
