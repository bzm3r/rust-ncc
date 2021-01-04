use crate::interactions::dat_utils::{sort_ixs, sym_sum};
use crate::utils::avro::Schematized;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};

#[derive(Clone, Default, Deserialize, Serialize, Schematize)]
pub struct SymCcDat<T: Copy + Default + Schematized> {
    pub num_cells: usize,
    sym_sum_nc: usize,
    dat: Vec<T>,
    undefined: T,
}

impl<T: Copy + Default + Schematized> SymCcDat<T> {
    pub fn new(num_cells: usize, undefined: T) -> SymCcDat<T> {
        let sym_sum_nc = sym_sum(num_cells);
        SymCcDat {
            num_cells,
            sym_sum_nc,
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
        let (small_ix, big_ix) = sort_ixs(ci, oci);
        sym_sum(self.num_cells) - sym_sum(self.num_cells - small_ix)
            + (big_ix - 1)
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
