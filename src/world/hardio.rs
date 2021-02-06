use crate::world::{Snapshot, WorldInfo};
use bincode::{deserialize_from, serialize_into};
use serde::Serialize;
use serde_cbor::ser::IoWrite;
use std::borrow::Borrow;
use std::fs::File;
use std::fs::OpenOptions;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{channel, Receiver, Sender};
use std::thread::JoinHandle;
use std::{io, thread};

#[allow(unused)]
#[derive(Clone, Copy)]
pub enum Format {
    Cbor,
    Bincode,
}

pub fn get_file_name(format: Format, name: &str) -> String {
    let ext = match format {
        Format::Cbor => "cbor".to_string(),
        Format::Bincode => "binc".to_string(),
    };
    if cfg!(features = "validate") {
        format!("validated_{}.{}", name, ext)
    } else {
        format!("{}.{}", name, ext)
    }
}

pub fn get_clean_file(
    out_dir: &PathBuf,
    name: &str,
    format: Format,
) -> io::Result<File> {
    let mut path = out_dir.clone();
    path.push(get_file_name(format, name));

    OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&path)
}

pub fn get_read_file(
    out_dir: &PathBuf,
    name: &str,
    format: Format,
) -> io::Result<File> {
    let mut path = out_dir.clone();
    path.push(get_file_name(format, name));

    OpenOptions::new().read(true).open(&path)
}

pub fn save_binc_to_cbor(binc_path: &PathBuf, cbor_path: &PathBuf) {
    let mut src =
        OpenOptions::new().read(true).open(binc_path).unwrap();
    let dst = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(cbor_path)
        .unwrap();
    let mut serializer =
        serde_cbor::Serializer::new(IoWrite::new(dst));
    let world_info: WorldInfo = deserialize_from(&mut src).unwrap();
    world_info.serialize(&mut serializer).unwrap();

    loop {
        let rd: bincode::Result<Vec<Snapshot>> =
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

pub struct AsyncWriter {
    pub output_dir: PathBuf,
    pub file_name: String,
    sender: Sender<Vec<Snapshot>>,
    buf: Vec<Snapshot>,
    max_capacity: usize,
    thread_handle: JoinHandle<()>,
    pub file_path: PathBuf,
}

impl AsyncWriter {
    pub fn new(
        output_dir: PathBuf,
        file_name: String,
        max_capacity: usize,
        truncate: bool,
        info: WorldInfo,
    ) -> AsyncWriter {
        let path = output_dir
            .join(get_file_name(Format::Bincode, &file_name));
        let (sender, receiver): (
            Sender<Vec<Snapshot>>,
            Receiver<Vec<Snapshot>>,
        ) = channel();
        let mut file = OpenOptions::new()
            .create(true)
            .write(true)
            .append(!truncate)
            .truncate(truncate)
            .open(&path)
            .unwrap();

        if truncate {
            serialize_into(&mut file, &info).unwrap();
        }

        let thread_handle = thread::spawn(move || {
            let mut f = file;
            let r = receiver;
            while let Ok(data_vec) = &r.recv() {
                serialize_into(&mut f, &data_vec).unwrap();
            }
        });

        AsyncWriter {
            output_dir,
            file_name,
            sender,
            buf: Vec::with_capacity(max_capacity),
            max_capacity,
            thread_handle,
            file_path: path,
        }
    }

    pub fn push(&mut self, data: Snapshot) {
        self.buf.push(data);
        if self.buf.len() == self.max_capacity {
            self.drain();
        }
    }

    pub fn drain(&mut self) {
        self.sender.send(self.buf.drain(..).collect()).unwrap();
    }

    pub fn finish(mut self, save_cbor: bool) {
        self.drain();
        let Self {
            sender,
            thread_handle,
            output_dir,
            file_name,
            file_path,
            ..
        } = self;
        drop(sender);
        thread_handle.join().unwrap();
        if save_cbor {
            let cbor_path = output_dir
                .join(get_file_name(Format::Cbor, &file_name));
            save_binc_to_cbor(&file_path, &cbor_path);
        }
    }
}

pub fn load(
    out_dir: &PathBuf,
    format: Format,
    name: &str,
) -> WorldInfo {
    let mut f = get_read_file(out_dir, name, format).unwrap();

    match format {
        Format::Cbor => serde_cbor::from_reader(&mut f).unwrap(),
        Format::Bincode => deserialize_from(&mut f).unwrap(),
    }
}

pub fn load_binc_from_path(file_path: &Path) -> File {
    if let Some(ext) = file_path.extension() {
        match ext.to_str().unwrap() {
            "binc" => {
                let f = OpenOptions::new()
                    .read(true)
                    .open(&file_path)
                    .unwrap();
                f
            }
            _ => panic!(
                "file path has unknown extension: {}",
                ext.to_str().unwrap()
            ),
        }
    } else {
        panic!(
            "no file extension in path: {}",
            file_path.to_str().unwrap()
        )
    }
}
