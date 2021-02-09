pub mod cell;
pub mod experiments;
pub mod hardio;
pub mod interactions;
pub mod math;
pub mod parameters;
pub mod utils;
pub mod world;

/// Number of vertices per model cell.
pub const NVERTS: usize = 16;
/// Default directory where simulation output will be placed.
pub const DEFAULT_OUTPUT_DIR: &str = "B:\\rust-ncc\\output";
