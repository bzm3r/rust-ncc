use crate::math::round;
use crate::NVERTS;

pub fn stringify_f64_arr(arr: &[f64; NVERTS], round_to: u32) -> String {
    arr.iter()
        .map(|&x| round(x, round_to).to_string())
        .collect::<Vec<String>>()
        .join(", ")
}
