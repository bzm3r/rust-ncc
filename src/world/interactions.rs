use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, Default, Deserialize, Schematize, Serialize)]
pub struct InteractionState {
    pub x_cil: [f32; NVERTS as usize],
    pub x_chemoas: [f32; NVERTS as usize],
    pub x_coas: [f32; NVERTS as usize],
    pub x_bdrys: [f32; NVERTS as usize],
}
