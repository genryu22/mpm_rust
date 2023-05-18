use std::sync::mpsc::{Receiver, Sender};

use mlsmpm::{Node, Particle};

use crate::Snapshot;

pub struct Executor {
    is_step_execution: bool,
    data_sender: Sender<Snapshot>,
    step_receiver: Receiver<usize>,
    target_step: Option<usize>,
}

impl Executor {
    pub fn new(
        data_sender: Sender<Snapshot>,
        step_receiver: Receiver<usize>,
        target_step: Option<usize>,
    ) -> Self {
        Self {
            is_step_execution: true,
            data_sender,
            step_receiver,
            target_step,
        }
    }

    pub fn start<T>(&mut self, mut step_executor: T)
    where
        T: StepExecutor,
    {
        let mut step = || {
            let data = step_executor.step();
            data
        };

        if let Some(target_step) = self.target_step {
            if target_step == 0 {
                return;
            }
        }

        loop {
            if self.is_step_execution {
                if let Ok(steps) = self.step_receiver.recv() {
                    if steps <= 0 {
                        self.is_step_execution = false;
                        continue;
                    }

                    if self.is_step_execution {
                        for _i in 0..steps {
                            let data = step();

                            if self.target_step.is_some() && data.steps >= self.target_step.unwrap()
                            {
                                self.data_sender.send(data).unwrap();
                                return;
                            } else if self.target_step.is_none() {
                                self.data_sender.send(data).unwrap();
                            }
                        }
                    }
                }
            } else {
                if let Ok(_) = self.step_receiver.try_recv() {
                    self.is_step_execution = true;
                    continue;
                }
                let data = step();
                if self.target_step.is_some() && data.steps >= self.target_step.unwrap() {
                    self.data_sender.send(data).unwrap();
                    return;
                } else if self.target_step.is_none() {
                    self.data_sender.send(data).unwrap();
                }
            }
        }
    }
}

pub trait StepExecutor {
    fn step(&mut self) -> Snapshot;
}
