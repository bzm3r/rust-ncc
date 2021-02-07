use bincode::{deserialize_from, serialize_into};
use rand::distributions::Uniform;
use rand::rngs::ThreadRng;
use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};
use std::borrow::Borrow;
use std::fmt::Display;
use std::fs::OpenOptions;
use std::path::Path;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::thread::JoinHandle;
use std::time::Duration;
use std::{fmt, io, thread};

pub const MAX_CAPACITY: usize = 2;
pub const DAT_LEN: usize = 2;
pub const SLEEP_TIME: Duration = Duration::from_millis(100);

#[derive(Clone, Copy, Deserialize, Serialize, PartialEq)]
pub struct Data {
    nums: [f32; DAT_LEN],
}

impl Data {
    fn random(rng: &mut ThreadRng, distrib: &Uniform<f32>) -> Data {
        Data {
            nums: {
                let mut r = [0.0_f32; DAT_LEN];
                rng.sample_iter(&distrib)
                    .take(DAT_LEN)
                    .enumerate()
                    .for_each(|(i, n)| r[i] = n);
                r
            },
        }
    }
}

impl Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}]",
            self.nums
                .iter()
                .map(|n| n.to_string())
                .collect::<Vec<String>>()
                .join(", ")
        )
    }
}

pub struct AsyncWriter {
    sender: Option<Sender<Vec<Data>>>,
    buf: Vec<Data>,
    max_capacity: usize,
    thread_handle: JoinHandle<()>,
}

impl AsyncWriter {
    pub fn new(
        path: &Path,
        max_capacity: usize,
        truncate: bool,
    ) -> AsyncWriter {
        let (sender, receiver): (
            Sender<Vec<Data>>,
            Receiver<Vec<Data>>,
        ) = channel();
        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .append(!truncate)
            .truncate(truncate)
            .open(path)
            .unwrap();

        let thread_handle = thread::spawn(move || {
            let mut f = file;
            let r = receiver;
            while let Ok(data_vec) = &r.recv() {
                println!(
                    "Received data: {}",
                    data_vec
                        .iter()
                        .map(|d| d.to_string())
                        .collect::<Vec<String>>()
                        .join(", ")
                );
                serialize_into(&mut f, &data_vec).unwrap();
            }
        });

        AsyncWriter {
            sender: Some(sender),
            buf: Vec::with_capacity(max_capacity),
            max_capacity,
            thread_handle,
        }
    }

    pub fn push(&mut self, data: Data) {
        self.buf.push(data);
        if self.buf.len() == self.max_capacity {
            self.drain();
        }
    }

    pub fn drain(&mut self) {
        println!(
            "Draining buffer: {}",
            self.buf
                .iter()
                .map(|d| d.to_string())
                .collect::<Vec<String>>()
                .join(", ")
        );
        if let Some(s) = self.sender.as_ref() {
            s.send(self.buf.drain(..).collect()).unwrap();
        }
    }

    pub fn finish(mut self) {
        self.drain();
        let Self {
            sender,
            thread_handle,
            ..
        } = self;
        drop(sender);
        thread_handle.join().unwrap();
    }
}

fn main() {
    let path = Path::new("play.dat");

    let mut rng = thread_rng();
    let uniform = Uniform::from(0.0_f32..1.0);
    let write_data: Vec<Data> =
        (0..25).map(|_| Data::random(&mut rng, &uniform)).collect();

    println!("Creating writer...");
    let mut writer = AsyncWriter::new(path, MAX_CAPACITY, true);
    println!("Beginning write...");
    for dat in write_data.iter() {
        println!("Pushing data: {}", dat);
        writer.push(*dat);
    }
    println!("Finishing writer...");
    writer.finish();
    println!("Writing complete.");

    let mut f = OpenOptions::new().read(true).open(path).unwrap();
    let mut read_data: Vec<Data> = vec![];
    loop {
        let rd: bincode::Result<Vec<Data>> = deserialize_from(&mut f);
        match rd {
            Ok(mut data) => read_data.append(&mut data),
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
    println!("(write_data, read_data)");
    for n in 0..write_data.len() {
        if n < read_data.len() {
            println!("({}, {})", write_data[n], read_data[n]);
        } else {
            println!("({}, EMPTY)", write_data[n]);
        }
    }
    println!("write_data len: {}", write_data.len());
    println!(
        "read_data matches write_data: {}",
        write_data.len() == read_data.len()
            && write_data
                .iter()
                .zip(read_data.iter())
                .all(|(wd, rd)| wd == rd)
    );
}
