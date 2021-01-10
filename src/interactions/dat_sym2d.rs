use crate::interactions::dat_utils::{sort_ixs, sym_sum};
use serde::{Deserialize, Serialize};

#[derive(Clone, Default, Deserialize, Serialize)]
pub struct SymCcDat<T: Copy + Default> {
    pub num_cells: usize,
    unique_elements: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy + Default> SymCcDat<T> {
    pub fn new(num_cells: usize, undefined: T) -> SymCcDat<T> {
        let sym_sum_nc = sym_sum(num_cells);
        SymCcDat {
            num_cells,
            unique_elements: sym_sum_nc,
            dat: vec![undefined; sym_sum_nc],
            undefined,
        }
    }

    /// Suppose we have an `n`x`n` symmetric matrix,
    /// whose diagonal elements are not of interest
    /// to us.
    ///
    /// The total number of elements in the matrix is `n^2`.
    /// There are `n` elements on the diagonal.
    /// Of the elements that remain, half are guaranteed
    /// to be non-unique, due to symmetry of the matrix.
    /// Therefore the total number of possibly unique
    /// elements is:
    /// `T = (n^2 - n)/2`
    ///
    /// Suppose we are required to find the `(i, j)` element
    /// of the matrix. First sort `(i, j)` into a tuple
    /// `(s, b)`, where `s = min(i, j)` and `b = max(i, j)`.
    /// Then, we try to find the element in the matrix at
    /// row `s`, column `b`, assuming that the elements are
    /// stored linearly, omitting those that are on the diagonal
    /// and guaranteed to be non-unique.
    ///
    /// Assuming we are on row `s`, how many elements have been
    /// stored thus far in the linear arrangement?
    ///
    /// On the `s`th row, we have
    /// `1` diagonal element, and `s` non-unique elements. Thus,
    /// relative to the elements which came before.
    fn calc_ix(&self, ci: usize, oci: usize) -> usize {
        // Suppose we have num_cells cells. Then, we would have an nxn
        // interaction matrix.
        //
        // Total number of possibly unique elements is num_cells*(num_cells - 1)/2,
        // because interactions in this case are symmetric.
        //
        // e.g. if nc = 4, then sym_sum_nc = 4*(4-1)/2 = 6
        // Suppose we get (ci, oci). We can assume that ci < oci,
        // because if it were not, we can sort ci/oci so that they are.
        // Note that because ci != oci (there is no self-interaction),
        // the inequality ci < oci is strict.
        let (small_ix, big_ix) = sort_ixs(ci, oci);

        // From now on, we will refer to (ci, oci) as (small_ix, big_ix).
        // We know that we have passed small_ix number of rows in the
        // matrix so far. Thus, the remaining number of possibly
        // unique elements in the matrix is sym_sum(num_cells - small_ix).

        // println!("num_cells: {}", self.num_cells);
        // println!("small_ix, big_ixL {}, {}", small_ix, big_ix);
        // println!("unique elements: {}", self.unique_elements);
        // println!(
        //     "remaining elements: {}",
        //     sym_sum(self.num_cells - small_ix)
        // );
        // println!(
        //     "delta: big_ix - small_ix - 1: {}",
        //     big_ix - small_ix - 1
        // );
        // println!(
        //     "calculated ix: {}",
        //     self.unique_elements - sym_sum(self.num_cells - small_ix)
        //         + (big_ix - small_ix - 1)
        // );
        self.unique_elements - sym_sum(self.num_cells - small_ix)
            + (big_ix - small_ix - 1)
    }

    pub fn set(&mut self, ci: usize, oci: usize, x: T) {
        if ci != oci {
            let ix = self.calc_ix(ci, oci);
            // println!("calculated ix: {}", ix);
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
