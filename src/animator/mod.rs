use crate::math::v2d::V2D;
use crate::world::{Cells, Snapshot};
use crate::NVERTS;
use cairo::{Context, Format, ImageSurface};
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

fn set_background(context: &Context) {
    context.set_source_rgb(1.0, 1.0, 1.0);
    context.paint();
}

fn draw_cell_poly(context: &Context, cell_poly: &[V2D; NVERTS]) {
    context.set_source_rgb(0.0, 0.0, 0.0);
    context.move_to(cell_poly[0].x as f64, cell_poly[0].y as f64);
    cell_poly[1..].iter().for_each(|v| {
        context.line_to(v.x as f64, v.y as f64);
    });
    context.close_path();
    context.set_line_width(2.0);
    context.stroke();
}

fn create_mp4(
    data: &DrawingData,
    width: i32,
    height: i32,
    framerate: i32,
    output_path: &Path,
) {
    //let frames = history.len();
    let surface = ImageSurface::create(Format::ARgb32, width, height)
        .expect("Couldn't create surface");
    let context = Context::new(&surface);
    let mut child = Command::new("ffmpeg")
        .args(&[
            "-f",
            "rawvideo",
            "-pix_fmt",
            "bgra",
            "-s",
            &format!("{}x{}", width, height),
            "-i",
            "-",
            "-pix_fmt",
            "yuv420p",
            "-r",
            &format!("{}", framerate),
            "-y",
            output_path.to_str().unwrap(),
        ])
        .stdin(Stdio::piped())
        .spawn()
        .expect("failed to execute process");
    {
        // limited borrow of stdin
        let child_stdin =
            child.stdin.as_mut().expect("Failed toget stdin");

        (0..data.num_frames).for_each(|frame| {
            set_background(&context);
            for cell_poly in data.get_cell_polys(frame) {
                draw_cell_poly(&context, cell_poly)
            }
            surface
                .with_data(|buf| {
                    child_stdin
                        .write_all(buf)
                        .expect("Failed to write bytes")
                })
                .expect("Failed to get_data");
        });
    }
    child.wait().expect("child process wasn't running");
}

pub struct DrawingData {
    num_cells: usize,
    num_frames: usize,
    //time_strings: Vec<String>,
    cell_polys: Vec<[V2D; NVERTS]>,
}

impl DrawingData {
    pub fn from_history(
        history: &[&Cells],
        px_w: i32,
        px_h: i32,
        px_per_micron: f32,
    ) -> DrawingData {
        let num_cells = history[0].cell_states.len();
        let mut cell_polys: Vec<[V2D; NVERTS]> = vec![];
        for cells in history.iter() {
            for cell in cells.cell_states.iter() {
                let mut transformed_vs = [V2D::default(); NVERTS];
                transformed_vs
                    .iter_mut()
                    .zip(cell.core.poly.iter())
                    .for_each(|(new_v, old_v)| {
                        *new_v =
                            old_v.scale(px_per_micron).translate(
                                px_w as f32 * 0.5,
                                px_h as f32 * 0.5,
                            );
                    });
                cell_polys.push(transformed_vs);
            }
        }
        DrawingData {
            num_cells,
            num_frames: history.len(),
            cell_polys,
        }
    }

    pub fn get_cell_polys(&self, frame: usize) -> &[[V2D; NVERTS]] {
        let start = frame * self.num_cells;
        let end = start + self.num_cells;
        &self.cell_polys[start..end]
    }
}

pub fn create_animation(history: &[Snapshot], output_path: &Path) {
    let cells_history =
        history.iter().map(|h| &h.cells).collect::<Vec<&Cells>>();
    let data =
        DrawingData::from_history(&cells_history, 1280, 720, 2.0);
    create_mp4(&data, 1000, 1000, 30, output_path);
}
