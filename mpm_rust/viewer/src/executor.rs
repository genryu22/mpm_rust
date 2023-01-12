use std::sync::mpsc::{Receiver, Sender};

use mlsmpm::Particle;

pub struct Executor {
    is_step_execution: bool,
    particles_sender: Sender<Vec<Particle>>,
    step_receiver: Receiver<usize>,
}

impl Executor {
    pub fn new(particles_sender: Sender<Vec<Particle>>, step_receiver: Receiver<usize>) -> Self {
        Self {
            is_step_execution: true,
            particles_sender,
            step_receiver,
        }
    }

    pub fn start<T>(&mut self, mut step_executor: T)
    where
        T: StepExecutor,
    {
        loop {
            if self.is_step_execution {
                if let Ok(steps) = self.step_receiver.recv() {
                    if steps <= 0 {
                        self.is_step_execution = false;
                        continue;
                    }

                    if self.is_step_execution {
                        for _i in 0..steps {
                            let particles = step_executor.step();
                            self.particles_sender.send(particles).unwrap();
                        }
                    }
                }
            } else {
                if let Ok(_) = self.step_receiver.try_recv() {
                    self.is_step_execution = true;
                    continue;
                }
                let particles = step_executor.step();
                self.particles_sender.send(particles).unwrap();
            }
        }
    }
}

pub trait StepExecutor {
    fn step(&mut self) -> Vec<Particle>;
}
