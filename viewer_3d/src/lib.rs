use std::{sync::mpsc, thread};

use executor::*;
use mlsmpm_3d::*;

mod executor;
mod window_bevy;
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

pub fn run_window_bevy(space_size: f64, settings: Settings, space: Space) {
    println!("{:?}", settings);

    let (step_sender, step_receiver) = mpsc::channel();
    let (data_sender, data_receiver) = mpsc::channel();

    step_sender.send(0).unwrap();

    thread::spawn(move || {
        let calc = Calculator::new(&settings, space);
        data_sender
            .send(Snapshot::new(
                calc.get_particles().to_vec(),
                calc.get_grid(),
                0,
            ))
            .unwrap();
        let step_executor = MLSMPMExecutor::new(calc);
        let mut executor = Executor::new(data_sender, step_receiver);
        executor.start(step_executor);
    });

    window_bevy::run(data_receiver, |snapshot| {});
}
