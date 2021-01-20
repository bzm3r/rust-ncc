use bincode::deserialize_from;
use crossbeam::channel::{bounded, Receiver, Sender};
use druid::im::Vector;
use rust_ncc::world::Snapshot;
use std::borrow::Borrow;
use std::fs::{File, OpenOptions};
use std::io;
use std::io::{Seek, SeekFrom};
use std::path::PathBuf;

pub struct AsyncReader {
    file_path: PathBuf,
    file: Option<File>,
    prev_offsets: Vec<u64>,
    rx_from_app: Receiver<Msg>,
    tx_to_app: Sender<Msg>,
}

pub enum Msg {
    FromReader(Payload),
    FromApp(Request),
}
pub enum Payload {
    Payload(Vector<Snapshot>),
    EndOfData,
    StartOfData,
    NoFile,
    LoadError(io::Error),
    LoadSuccess,
}

pub enum Request {
    FetchNext,
    FetchPrev,
    LoadFile(PathBuf),
}

pub enum DeserializeError {
    EndOfFile,
    Other(bincode::Error),
}

impl AsyncReader {
    pub fn new(tx_to_app: Sender<Msg>) -> (Sender<Msg>, AsyncReader) {
        let (tx_to_reader, rx_from_app) = bounded(1);

        (
            tx_to_reader,
            AsyncReader {
                file_path: Default::default(),
                file: None,
                prev_offsets: vec![],
                rx_from_app,
                tx_to_app,
            },
        )
    }

    fn fetch_next(&mut self) {
        if let Some(f) = self.file.as_mut() {
            self.prev_offsets
                .push(f.seek(SeekFrom::Current(0)).unwrap());
            let snaps: bincode::Result<Vec<Snapshot>> =
                deserialize_from(f);
            match snaps {
                Ok(snaps) => {
                    self.tx_to_app
                        .send(Msg::FromReader(Payload::Payload(
                            Vector::from(&snaps),
                        )))
                        .unwrap();
                }
                Err(err) => {
                    if let bincode::ErrorKind::Io(std_err) =
                        err.borrow()
                    {
                        if let io::ErrorKind::UnexpectedEof =
                            std_err.kind()
                        {
                            self.tx_to_app
                                .send(Msg::FromReader(
                                    Payload::EndOfData,
                                ))
                                .unwrap();
                        }
                    }
                    panic!(err);
                }
            }
        } else {
            self.tx_to_app
                .send(Msg::FromReader(Payload::NoFile))
                .unwrap();
        }
    }

    fn fetch_prev(&mut self) {
        if let Some(f) = self.file.as_mut() {
            match self.prev_offsets.pop() {
                None => self
                    .tx_to_app
                    .send(Msg::FromReader(Payload::StartOfData))
                    .unwrap(),
                Some(offset) => {
                    f.seek(SeekFrom::Start(offset)).unwrap();
                    let snaps: bincode::Result<Vec<Snapshot>> =
                        deserialize_from(f);
                    match snaps {
                        Ok(snaps) => {
                            self.tx_to_app
                                .send(Msg::FromReader(
                                    Payload::Payload(Vector::from(
                                        &snaps,
                                    )),
                                ))
                                .unwrap();
                        }
                        Err(err) => {
                            if let bincode::ErrorKind::Io(std_err) =
                                err.borrow()
                            {
                                if let io::ErrorKind::UnexpectedEof =
                                    std_err.kind()
                                {
                                    self.tx_to_app
                                        .send(Msg::FromReader(
                                            Payload::StartOfData,
                                        ))
                                        .unwrap();
                                }
                            }
                            panic!(err);
                        }
                    }
                }
            }
        } else {
            self.tx_to_app
                .send(Msg::FromReader(Payload::NoFile))
                .unwrap();
        }
    }

    fn load_file(&mut self, path: PathBuf) {
        match OpenOptions::new().read(true).open(&path) {
            Ok(f) => {
                self.file = Some(f);
                self.file_path = path;
                self.prev_offsets.truncate(0);
                self.tx_to_app
                    .send(Msg::FromReader(Payload::LoadSuccess))
                    .unwrap();
            }
            Err(err) => {
                self.tx_to_app
                    .send(Msg::FromReader(Payload::LoadError(err)))
                    .unwrap();
            }
        }
    }

    pub fn work_loop(&mut self) {
        loop {
            match self.rx_from_app.recv() {
                Ok(msg) => match msg {
                    Msg::FromApp(Request::FetchNext) => {
                        self.fetch_next();
                    }
                    Msg::FromApp(Request::FetchPrev) => {
                        self.fetch_prev();
                    }
                    Msg::FromApp(Request::LoadFile(path)) => {
                        self.load_file(path)
                    }
                    Msg::FromReader(_) => {
                        panic!("Received message from reader, but this is the reader.");
                    }
                },
                Err(_recv_err) => {
                    break;
                }
            }
        }
    }
}
