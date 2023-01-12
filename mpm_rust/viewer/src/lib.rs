use std::{
    sync::{
        mpsc::{self, Receiver, Sender},
        Arc, Mutex,
    },
    thread,
    time::Duration,
};

pub mod executor;
pub mod window;

use mlsmpm::*;

use crate::executor::Executor;

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
    fn step(&mut self) -> Vec<Particle> {
        self.calculator.update();
        self.steps += 1;
        self.calculator.get_particles().to_vec()
    }
}

pub fn run_window(space_size: f64, settings: Settings, space: Space) {
    println!("{:?}", settings);

    let (step_sender, step_receiver) = mpsc::channel();
    let (particles_sender, particles_receiver) = mpsc::channel();

    thread::spawn(move || {
        let calc = Calculator::new(&settings, space);
        particles_sender
            .send(calc.get_particles().to_vec())
            .unwrap();
        let step_executor = MLSMPMExecutor::new(calc);
        let mut executor = Executor::new(particles_sender, step_receiver);
        executor.start(step_executor);
    });

    let mut window = window::ParticleWindow::new(space_size, particles_receiver, step_sender);
    window.run();
}

pub fn run_dambreak_window() {
    let settings = Settings {
        dt: 1e-3,
        gravity: -100.,
        dynamic_viscosity: 1e-6,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 100,
        c: 1e2,
        eos_power: 4.,
    };

    let space = Space::new_for_dambreak(&settings);
    run_window(settings.space_width, settings, space);
}

pub fn run_poiseuille_window() {
    let settings = Settings {
        dt: 0.01,
        gravity: -1e-2,
        dynamic_viscosity: 1e-2,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 200,
        c: 0.,
        eos_power: 0.,
    };

    let v_time_steps = (1. / settings.dynamic_viscosity / settings.dt).ceil() as u32;
    println!("粘性時間: L^2/mu = {} steps", v_time_steps);
    println!(
        "平均流速の理論解: U_mean = g*L*L/(nu*12) = {}",
        settings.gravity * 1. * 1. / ((settings.dynamic_viscosity / 1.) * 12.)
    );
    println!(
        "最大流速の理論解: U_max = 1.5*U_mean = {}",
        1.5 * settings.gravity * 1. * 1. / ((settings.dynamic_viscosity / 1.) * 12.)
    );

    let space = Space::new_for_poiseuille(&settings);
    run_window(1., settings, space);
}
