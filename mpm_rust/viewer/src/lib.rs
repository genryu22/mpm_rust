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

pub fn run_dambreak_window_wgpu() {
    let settings = Settings {
        dt: 1e-3,
        gravity: -100.,
        dynamic_viscosity: 1e-3,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 100,
        rho_0: 1000.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: true,
        vx_zero: false,
    };

    let space = Space::new_for_dambreak(&settings);
    run_window_wgpu(settings.space_width, settings, space);
}

pub fn run_window_bevy(space_size: f64, settings: Settings, space: Space) {
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

    window_bevy::run(data_receiver);
}

pub fn run_dambreak_window_bevy() {
    let settings = Settings {
        dt: 1e-3,
        gravity: -100.,
        dynamic_viscosity: 1e-3,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 100,
        rho_0: 1000.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: true,
        vx_zero: false,
    };

    let space = Space::new_for_dambreak(&settings);
    run_window_bevy(settings.space_width, settings, space);
}

pub fn run_window(space_size: f64, settings: Settings, space: Space) {
    println!("{:?}", settings);

    let (step_sender, step_receiver) = mpsc::channel();
    let (data_sender, data_receiver) = mpsc::channel();

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

    let mut window = window::ParticleWindow::new(space_size, data_receiver, step_sender);
    window.run(&|snapshot| write_to_files(snapshot).unwrap());
}

pub fn run_dambreak_window() {
    let settings = Settings {
        dt: 1e-3,
        gravity: -100.,
        dynamic_viscosity: 1e-3,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 100,
        rho_0: 1000.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: true,
        vx_zero: false,
    };

    let space = Space::new_for_dambreak(&settings);
    run_window(settings.space_width, settings, space);
}

pub fn run_dambreak_experiment_window() {
    let settings = Settings {
        dt: 1e-4,
        gravity: -10.,
        dynamic_viscosity: 1e-3,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 500,
        rho_0: 1000.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: true,
        vx_zero: false,
    };

    let space = Space::new_for_dambreak_experiment(&settings);
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
        rho_0: 1.,
        c: 0.,
        eos_power: 0.,
        boundary_mirror: true,
        vx_zero: true,
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

pub fn run_taylorgreen_window() {
    let settings = Settings {
        dt: 4e-4,
        gravity: 0.,
        dynamic_viscosity: 1e-3,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 200,
        rho_0: 1.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: false,
        vx_zero: false,
    };

    let space = Space::new_for_taylor_green(&settings);
    run_window(7., settings, space);
}
