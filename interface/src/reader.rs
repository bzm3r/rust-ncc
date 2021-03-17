use bincode::deserialize_from;
use crossbeam::channel::{bounded, Sender};
use druid::im::Vector;
use druid::{ExtEventSink, Target};
use rust_ncc::hardio::WorldSnapshot;
use std::borrow::Borrow;
use std::fs::{File, OpenOptions};
use std::io::{Seek, SeekFrom};
use std::path::PathBuf;
use std::{io, thread};

pub struct AsyncReader {
    file_path: PathBuf,
    file: Option<File>,
    prev_offsets: Vec<u64>,
}

pub struct AppContact {
    request: FetchRequest,
    event_sink: ExtEventSink,
}

impl Default for FetchRequest {
    fn default() -> Self {
        FetchRequest::NoOp
    }
}

#[derive(Clone)]
pub enum FetchResult {
    Payload(Vector<WorldSnapshot>),
    FileBoundsHit,
    NoFile,
    IOError(io::ErrorKind),
    BincodeError(String),
    LoadSuccess,
    NoOp,
}

impl From<bincode::Error> for FetchResult {
    fn from(err: bincode::Error) -> Self {
        if let bincode::ErrorKind::Io(std_err) = err.borrow() {
            if let io::ErrorKind::UnexpectedEof = std_err.kind() {
                FetchResult::FileBoundsHit
            } else {
                FetchResult::IOError(std_err.kind())
            }
        } else {
            FetchResult::BincodeError(format!("{}", err))
        }
    }
}

#[derive(Clone)]
pub enum FetchRequest {
    NoOp,
    FetchNext,
    FetchPrev,
    LoadFile(PathBuf),
}

pub enum DeserializeError {
    EndOfFile,
    Other(bincode::Error),
}

impl AsyncReader {
    pub fn spawn_for_path(path: &PathBuf) -> Sender<AppContact> {
        let (tx_to_reader, rx_from_app) = bounded::<AppContact>(1);
        let mut reader = AsyncReader {
            file_path: path.clone(),
            file: None,
            prev_offsets: vec![],
        };

        thread::spawn(move || loop {
            match rx_from_app.recv() {
                Ok(fetch_request) => reader.handle(fetch_request),
                Err(_) => break,
            }
        });

        tx_to_reader
    }

    fn store_file_offset(&mut self, f: &mut File) {
        self.prev_offsets
            .push(f.seek(SeekFrom::Current(0)).unwrap());
    }

    fn fetch_next(&mut self) -> FetchResult {
        if let Some(f) = self.file.as_mut() {
            self.store_file_offset(f);
            let snaps: bincode::Result<Vec<WorldSnapshot>> =
                deserialize_from(f);
            match snaps {
                Ok(snaps) => {
                    FetchResult::Payload(Vector::from(&snaps))
                }
                Err(err) => FetchResult::from(err),
            }
        } else {
            FetchResult::NoFile
        }
    }

    fn fetch_prev(&mut self) -> FetchResult {
        if let Some(f) = self.file.as_mut() {
            match self.prev_offsets.pop() {
                None => FetchResult::FileBoundsHit,
                Some(offset) => {
                    f.seek(SeekFrom::Start(offset)).unwrap();
                    let snaps: bincode::Result<Vec<WorldSnapshot>> =
                        deserialize_from(f);
                    match snaps {
                        Ok(snaps) => {
                            FetchResult::Payload(Vector::from(&snaps))
                        }
                        Err(err) => FetchResult::from(err),
                    }
                }
            }
        } else {
            FetchResult::NoFile
        }
    }

    fn load_file(&mut self, path: PathBuf) -> FetchResult {
        match OpenOptions::new().read(true).open(&path) {
            Ok(f) => {
                self.file = Some(f);
                self.file_path = path;
                self.prev_offsets.truncate(0);
                FetchResult::LoadSuccess
            }
            Err(err) => FetchResult::IOError(err.kind()),
        }
    }

    pub fn handle(&mut self, fetch_request: AppContact) {
        let AppContact {
            request,
            event_sink,
        } = fetch_request;
        event_sink.submit_command(
            super::delegate::FETCH_FULFILL,
            match request {
                FetchRequest::FetchPrev => self.fetch_prev(),
                FetchRequest::FetchNext => self.fetch_next(),
                FetchRequest::LoadFile(path) => self.load_file(path),
                FetchRequest::NoOp => FetchResult::NoOp,
            },
            Target::Global,
        );
    }
}
