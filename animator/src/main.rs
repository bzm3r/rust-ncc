use cairo::{Context, Format, ImageSurface};
use std::f64::consts::PI;
use std::io::Write;
use std::process::{Command, Stdio};

// fn draw(surface: &ImageSurface, f: i32) {
//     let cr = Context::new(surface);
//     let width = surface.get_width() as f64;
//     let height = surface.get_height() as f64;
//     let dx = f as f64 * 1.0;
//
//     cr.set_source_rgb(1.0, 1.0, 1.0);
//     cr.paint();
//
//     let cx = width / 2.0;
//     let cy = height / 2.0;
//     let r = cy;
//     let cstart = -0.5 * PI;
//     let cend = cstart + 2.0 * PI * ((f + 1) as f64) / 300.0;
//     cr.move_to(cx, cy);
//     cr.line_to(cx, 0.0);
//     cr.arc(cx, cy, r, cstart, cend);
//     cr.line_to(cx, cy);
//     cr.set_source_rgba(0.0, 0.5, 0.0, 0.2);
//     cr.fill();
//
//     cr.select_font_face(
//         "sans-serif",
//         cairo::FontSlant::Normal,
//         cairo::FontWeight::Normal,
//     );
//     cr.set_font_size(70.0);
//     cr.move_to(600.0 - dx, 100.0);
//     cr.set_source_rgb(0.0, 0.0, 1.0);
//     cr.show_text("Hello, world! 1234567890");
//     cr.fill();
// }
//
// fn make_movie(width: i32, height: i32, framerate: i32, frames: i32) {
//     let mut surface =
//         ImageSurface::create(Format::ARgb32, width, height).expect("Couldn't create surface");
//     let mut child = Command::new("ffmpeg")
//         .args(&[
//             "-f",
//             "rawvideo",
//             "-pix_fmt",
//             "bgra",
//             "-s",
//             &format!("{}x{}", width, height),
//             "-i",
//             "-",
//             "-pix_fmt",
//             "yuv420p",
//             "-r",
//             &format!("{}", framerate),
//             "-y",
//             "out.mp4",
//         ])
//         .stdin(Stdio::piped())
//         .spawn()
//         .expect("failed to execute process");
//     {
//         // limited borrow of stdin
//         let child_stdin = child.stdin.as_mut().expect("failed to get stdin");
//
//         (0..frames).for_each(|f| {
//             draw(&surface, f);
//             let d = surface.get_data().expect("Failed to get_data");
//             child_stdin.write_all(&d).expect("Failed to write to file");
//         });
//     }
//     child.wait().expect("child process wasn't running");
// }
//
fn main() {
    let output_path = std::path::Path::new("../output");
    make_movie(1280, 720, 30, 300);
}
