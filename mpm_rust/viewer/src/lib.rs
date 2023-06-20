use std::{
    sync::{
        mpsc::{self, Receiver, Sender},
        Arc, Mutex,
    },
    thread,
    time::Duration,
    vec,
};

pub mod executor;
mod file;
pub mod window;
pub mod window_bevy;
pub mod window_wgpu;

use mlsmpm::*;

use crate::{executor::Executor, file::write_to_files};

pub struct Snapshot {
    particles: Vec<Particle>,
    grid: Vec<Node>,
    steps: usize,
}

impl Snapshot {
    pub fn new(particles: Vec<Particle>, grid: Vec<Node>, steps: usize) -> Self {
        Self {
            particles,
            grid,
            steps,
        }
    }

    pub fn empty() -> Self {
        Self {
            particles: vec![],
            grid: vec![],
            steps: 0,
        }
    }
}

pub struct MLSMPMExecutor<'a> {
    steps: usize,
    calculator: Calculator<'a>,
}

impl<'a> MLSMPMExecutor<'a> {
    pub fn new(calculator: Calculator<'a>) -> Self {
        Self {
            steps: 0,
            calculator,
        }
    }
}

impl<'a> executor::StepExecutor for MLSMPMExecutor<'a> {
    fn step(&mut self) -> Snapshot {
        self.calculator.update();
        self.steps += 1;

        Snapshot::new(
            self.calculator.get_particles().to_vec(),
            self.calculator.get_grid().to_vec(),
            self.steps,
        )
    }
}

pub fn run_window_wgpu(space_size: f64, settings: Settings, space: Space) {
    println!("{:?}", settings);

    let count = space.get_particle_count() as u32;

    let (step_sender, step_receiver) = mpsc::channel();
    let (data_sender, data_receiver) = mpsc::channel();

    step_sender.send(0).unwrap();

    thread::spawn(move || {
        let calc = Calculator::new(&settings, space);
        data_sender
            .send(Snapshot::new(
                calc.get_particles().to_vec(),
                calc.get_grid().to_vec(),
                0,
            ))
            .unwrap();
        let step_executor = MLSMPMExecutor::new(calc);
        let mut executor = Executor::new(data_sender, step_receiver);
        executor.start(step_executor);
    });

    pollster::block_on(window_wgpu::run(count, data_receiver));
}

pub fn run_window_bevy(space_size: f64, settings: Settings, space: Space, window_size: f32) {
    println!("{:?}", settings);

    let count = space.get_particle_count() as u32;

    let (step_sender, step_receiver) = mpsc::channel();
    let (data_sender, data_receiver) = mpsc::channel();

    step_sender.send(0).unwrap();

    thread::spawn(move || {
        let calc = Calculator::new(&settings, space);
        data_sender
            .send(Snapshot::new(
                calc.get_particles().to_vec(),
                calc.get_grid().to_vec(),
                0,
            ))
            .unwrap();
        let step_executor = MLSMPMExecutor::new(calc);
        let mut executor = Executor::new(data_sender, step_receiver);
        executor.start(step_executor);
    });

    window_bevy::run(
        data_receiver,
        |snapshot| {
            if false {
                write_to_files(snapshot).unwrap()
            }
        },
        window_size,
    );
}
