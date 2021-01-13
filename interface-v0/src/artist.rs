use crate::polygon::Polygon;
use druid::{Data, Env, PaintCtx, Point, TextLayout, UpdateCtx};
use rust_ncc::parameters::quantity::Quantity;
use rust_ncc::world::History;
use std::sync::Arc;

#[derive(Clone, Copy, Data)]
pub struct Scales {
    pub rgtp: f32,
    pub force: f32,
    pub crl: f32,
}

#[derive(Clone, Data)]
pub struct Artist {
    time_in_secs: f32,
    polys: Arc<Vec<Polygon>>,
    scales: Scales,
}

impl Artist {
    pub fn update(
        &mut self,
        tstep: u32,
        h: &History,
        scales: Scales,
    ) -> Artist {
        let snapshot = &h.snapshots[tstep as usize];
        let time_in_secs =
            h.char_quants.t.mul_number(tstep as f32).number();
        let mut polys = Arc::new(
            snapshot
                .cells
                .states
                .iter()
                .zip(snapshot.cells.interactions.iter())
                .map(|(cell, interactions)| {
                    Polygon::new(
                        cell,
                        interactions,
                        &h.cell_params[0],
                        scales,
                    )
                })
                .collect::<Vec<Polygon>>(),
        );

        Artist {
            time_in_secs,
            polys,
            scales,
        }
    }

    pub fn paint(ctx: &mut PaintCtx, data: &Artist, env: &Env) {
        for poly in data.polys.iter() {
            Polygon::paint(ctx, poly, env);
        }
        let mut time = TextLayout::default();
        time.set_text(format!(
            "time = {} min.",
            (data.time_in_secs / 60.0) as u32
        ));
        time.draw(ctx, Point::ORIGIN);
    }
}

impl Default for Artist {
    fn default() -> Self {
        Artist {
            time_in_secs: 0.0,
            polys: Default::default(),
            scales: Scales {
                rgtp: 0.0,
                force: 0.0,
                crl: 0.0,
            },
        }
    }
}

//
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
