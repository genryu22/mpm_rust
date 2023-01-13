use std::sync::mpsc::{Receiver, Sender};

use mlsmpm::{Node, Particle};

use crate::Snapshot;

pub struct Executor {
    is_step_execution: bool,
    data_sender: Sender<Snapshot>,
    step_receiver: Receiver<usize>,
}

impl Executor {
    pub fn new(data_sender: Sender<Snapshot>, step_receiver: Receiver<usize>) -> Self {
        Self {
            is_step_execution: true,
            data_sender,
            step_receiver,
        }
    }

    pub fn start<T>(&mut self, mut step_executor: T)
    where
        T: StepExecutor,
    {
        let mut step = || {
            let data = step_executor.step();
            self.data_sender.send(data).unwrap();
        };

        loop {
            if self.is_step_execution {
                if let Ok(steps) = self.step_receiver.recv() {
                    if steps <= 0 {
                        self.is_step_execution = false;
                        continue;
                    }

                    if self.is_step_execution {
                        for _i in 0..steps {
                            step();
                        }
                    }
                }
            } else {
                if let Ok(_) = self.step_receiver.try_recv() {
                    self.is_step_execution = true;
                    continue;
                }
                step();
            }
        }
    }
}

pub trait StepExecutor {
    fn step(&mut self) -> Snapshot;
}
