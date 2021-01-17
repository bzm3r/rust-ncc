use bincode::{deserialize_from, serialize_into};
use rand::distributions::Uniform;
use rand::rngs::ThreadRng;
use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::time::Duration;
use std::{fmt, thread};

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
    pub receiver: Receiver<Vec<Data>>,
    pub file: File,
}

impl AsyncWriter {
    pub fn new(
        path: &Path,
        recv: Receiver<Vec<Data>>,
    ) -> AsyncWriter {
        AsyncWriter {
            receiver: recv,
            file: OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(path)
                .unwrap(),
        }
    }

    pub fn write(mut self) {
        thread::sleep(SLEEP_TIME);
        while let Ok(data_vec) = self.receiver.recv() {
            println!(
                "received: {}",
                data_vec
                    .iter()
                    .map(|d| d.to_string())
                    .collect::<Vec<String>>()
                    .join(", ")
            );
            serialize_into(&mut self.file, &data_vec).unwrap();
        }
    }
}

fn main() {
    let path = Path::new("play.dat");

    let mut rng = thread_rng();
    let uniform = Uniform::from(0.0_f32..1.0);
    let write_data: Vec<Data> =
        (0..5).map(|_| Data::random(&mut rng, &uniform)).collect();
    let mut send_data = vec![];

    let (sender, receiver): (Sender<Vec<Data>>, Receiver<Vec<Data>>) =
        channel();
    let async_writer = AsyncWriter::new(path, receiver);

    let writer = thread::spawn(move || async_writer.write());

    for n in 0..write_data.len() {
        send_data.push(write_data[n]);
        if send_data.len() == MAX_CAPACITY
            || n == (write_data.len() - 1)
        {
            let sd: Vec<Data> = send_data.drain(..).collect();
            println!(
                "sending: {}",
                sd.iter()
                    .map(|d| d.to_string())
                    .collect::<Vec<String>>()
                    .join(", ")
            );
            sender.send(sd).unwrap();
        }
    }
    drop(sender);

    writer.join().unwrap();
    let mut f = OpenOptions::new().read(true).open(path).unwrap();
    let read_data: Vec<Data> = deserialize_from(&mut f).unwrap();
    println!("(write_data, read_data)");
    for n in 0..write_data.len() {
        if n < read_data.len() {
            println!("({}, {})", write_data[n], read_data[n]);
        } else {
            println!("({}, EMPTY)", write_data[n]);
        }
    }
    println!(
        "read_data matches write_data: {}",
        write_data.len() == read_data.len()
            && write_data
                .iter()
                .zip(read_data.iter())
                .all(|(wd, rd)| wd == rd)
    );
}
