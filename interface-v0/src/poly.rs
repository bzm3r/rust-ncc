// use crate::canvas::Scales;
use druid::Color;
// use rust_ncc::cell::Cell;
// use rust_ncc::interactions::Interactions;
// use rust_ncc::math::geometry::Poly;
// use rust_ncc::parameters::Parameters;
// use rust_ncc::utils::circ_ix_plus;
// use rust_ncc::NVERTS;
use std::cmp::Ordering;

//
// fn gen_rgtp_edge_colors(
//     cell: &Cell,
//     scale: f32,
//     parameters: &Parameters,
// ) -> [LinearGradient; NVERTS] {
//     (0..NVERTS).for_each(|ix| {
//         let next_ix = circ_ix_plus(ix, NVERTS);
//         rgtp_edge_colors[ix] = LinearGradient::new(
//             UnitPoint::new(0.0, 0.5),
//             UnitPoint::new(1.0, 0.5),
//             vec![
//                 map_to_color(
//                     cell.core.rho_acts[ix],
//                     cell.core.rac_acts[ix],
//                     Some(relative_rgtp_activity[ix].to_f32()),
//                     Color::RED,
//                     Color::BLUE,
//                     scale,
//                 ),
//                 map_to_color(
//                     cell.core.rho_acts[next_ix],
//                     cell.core.rac_acts[next_ix],
//                     Some(relative_rgtp_activity[next_ix].to_f32()),
//                     Color::RED,
//                     Color::BLUE,
//                     scale,
//                 ),
//             ]
//             .into(),
//         );
//     });
//     rgtp_edge_colors
// }
//
// fn gen_crl_edge_colors(
//     interactions: &Interactions,
//     scale: f32,
// ) -> [LinearGradient; NVERTS] {
//     let mut crl_edge_colors = [default_lin_grad(); NVERTS];
//     (0..NVERTS).for_each(|ix| {
//         let next_ix = circ_ix_plus(ix, NVERTS);
//         crl_edge_colors[ix] = LinearGradient::new(
//             UnitPoint::new(0.0, 0.5),
//             UnitPoint::new(1.0, 0.5),
//             vec![
//                 map_to_color(
//                     interactions.x_cils[ix],
//                     interactions.x_cals[ix],
//                     None,
//                     Color::YELLOW,
//                     Color::AQUA,
//                     scale,
//                 ),
//                 map_to_color(
//                     interactions.x_cals[next_ix],
//                     interactions.x_cils[next_ix],
//                     None,
//                     Color::YELLOW,
//                     Color::AQUA,
//                     scale,
//                 ),
//             ]
//             .into(),
//         )
//     });
//     crl_edge_colors
// }
//
// fn gen_coa_edge_colors(
//     interactions: &Interactions,
//     scale: f32,
// ) -> [LinearGradient; NVERTS] {
//     let mut coa_edge_colors = [default_lin_grad(); NVERTS];
//     (0..NVERTS).for_each(|ix| {
//         let next_ix = circ_ix_plus(ix, NVERTS);
//         coa_edge_colors[ix] = LinearGradient::new(
//             UnitPoint::new(0.0, 0.5),
//             UnitPoint::new(1.0, 0.5),
//             vec![
//                 map_to_color(
//                     interactions.x_coas[ix],
//                     0.0,
//                     None,
//                     Color::BLACK,
//                     Color::FUCHSIA,
//                     scale,
//                 ),
//                 map_to_color(
//                     interactions.x_coas[next_ix],
//                     0.0,
//                     None,
//                     Color::BLACK,
//                     Color::FUCHSIA,
//                     scale,
//                 ),
//             ]
//             .into(),
//         )
//     });
//     coa_edge_colors
// }
//
