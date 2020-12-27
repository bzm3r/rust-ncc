use crate::NVERTS;

/// Given a pair of cell indices `(ci, oci)` (`ci` for "cell index",
/// and `oci` for "other cell index), return them in ordered format
/// `(s, b)` (`s` for smaller, and `b` for bigger), where `s = min(ci,
/// oci)` and `b = max(ci, oci)`.
pub fn sort_ixs(ci: usize, oci: usize) -> (usize, usize) {
    if oci < ci {
        (oci, ci)
    } else {
        (ci, oci)
    }
}

/// In a 4D data matrix, indexed by `(ci, vi, oci, ovi)`, where  
pub fn fix_oci(ci: usize, oci: usize) -> usize {
    if ci < oci {
        oci - 1
    } else {
        oci
    }
}

/// Given a cell-cell vertex-vertex index of the form `(ci, vi, oci,
/// ovi)` , where `vi` (`vi` for "vertex index") is the index of a
/// vertex on cell `ci` (`ci` for "cell index"), and `ovi` (`ovi` for
/// "other vertex index") is the index of a vertex on `oci` (`oci` for
/// "other cell index"), return `(sc, sv, bc, bv)`, where: `sc =
/// min(ci, oci)`, `bc = min(ci, oci)`, `sv` is the `vi` if `sc ==
/// ci`, or `sv == ovi`, if `sc == oci`. Similarly for `bv`.
pub fn sort_ixs4d(
    ci: usize,
    vi: usize,
    oci: usize,
    ovi: usize,
) -> (usize, usize, usize, usize) {
    if oci < ci {
        (oci, ovi, ci, vi)
    } else {
        (ci, vi, oci, ovi)
    }
}

// Debug mode utility function to confirm that
// #[allow(unused)]
// #[cfg(feature = "custom_debug")]
// pub fn check_indices4d(
//     num_cells: usize,
//     ci: usize,
//     vi: usize,
//     oci: usize,
//     ovi: usize,
// ) {
//     let oci = fix_oci(ci, oci);
//
//     if ci > num_cells {
//         panic!("{} cells tracked, received ci: {}", num_cells, ci);
//     }
//
//     if vi > NVERTS as usize {
//         panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
//     }
//
//     if oci > num_cells - 1 {
//         panic!(
//             "{} cells tracked, received oci: {} (max: {})",
//             num_cells,
//             oci,
//             num_cells - 1
//         );
//     }
//
//     if ovi > NVERTS as usize {
//         panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
//     }
// }

/// Given an `n`x`n` matrix, calculate the number of elements left are
/// removing the `n` diagonal elements, and de-duplication of the
/// remainder (`/2)`: `R = (n^2 - n) / 2 = n * (n - 1) / 2`
pub fn sym_sum(n: usize) -> usize {
    n * (n - 1) / 2
}

// #[cfg(feature = "custom_debug")]
// pub fn check_sym_indices4d(
//     num_cells: usize,
//     ci: usize,
//     vi: usize,
//     oci: usize,
//     ovi: usize,
// ) {
//     let (ci, vi, oci, ovi) = sort_ixs4d(ci, vi, oci, ovi);
//
//     if ci > num_cells {
//         panic!("{} cells tracked, received ci: {}", num_cells, ci);
//     }
//
//     if vi > NVERTS as usize {
//         panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
//     }
//
//     if oci > num_cells {
//         panic!(
//             "{} cells tracked, received symmetric oci: {}",
//             num_cells, oci
//         );
//     }
//
//     if ovi > NVERTS as usize {
//         panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
//     }
// }
