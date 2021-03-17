use bincode::deserialize_from;
use serde::{Deserialize, Serialize};
use serde_cbor::ser::IoWrite;
use simulator::world::{WorldInfo, WorldSnapshot};
use std::borrow::Borrow;
use std::fs::OpenOptions;
use std::io;

const BINC_PATH: &str = "../output/history_n_cells.binc";
const CBOR_PATH: &str = "../output/history_n_cells.cbor";

#[derive(Deserialize, Serialize, Debug)]
pub enum CborData {
    WorldInfo(WorldInfo),
    Snapshot(WorldSnapshot),
}

fn cbor_write() {
    let mut src =
        OpenOptions::new().read(true).open(BINC_PATH).unwrap();
    let dst = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(CBOR_PATH)
        .unwrap();
    let mut serializer =
        serde_cbor::Serializer::new(IoWrite::new(dst));
    let world_info: WorldInfo = deserialize_from(&mut src).unwrap();
    world_info.serialize(&mut serializer).unwrap();

    loop {
        let rd: bincode::Result<Vec<WorldSnapshot>> =
            deserialize_from(&mut src);
        match rd {
            Ok(snaps) => {
                snaps.serialize(&mut serializer).unwrap();
            }
            Err(err) => {
                if let bincode::ErrorKind::Io(std_err) = err.borrow()
                {
                    if let io::ErrorKind::UnexpectedEof =
                        std_err.kind()
                    {
                        break;
                    }
                }
                panic!("{}", err);
            }
        }
    }
}

fn cbor_read() {
    let mut src =
        OpenOptions::new().read(true).open(CBOR_PATH).unwrap();
    let mut deserializer =
        serde_cbor::Deserializer::from_reader(&mut src);
    let _world_info: WorldInfo =
        serde::Deserialize::deserialize(&mut deserializer).unwrap();

    let mut snapshots: Vec<WorldSnapshot> = vec![];
    loop {
        let rd: serde_cbor::error::Result<Vec<WorldSnapshot>> =
            serde::Deserialize::deserialize(&mut deserializer);
        match rd {
            Ok(mut snaps) => snapshots.append(&mut snaps),
            Err(err) => {
                if err.is_eof() {
                    break;
                } else {
                    panic!("{}", err);
                }
            }
        }
    }
    println!("num_snapshots: {}", snapshots.len());
}

fn main() {
    cbor_write();
    cbor_read();

    // println!("num_snapshots: {}", snapshots.len());
}
