use std::{
    sync::{
        mpsc::{Receiver, Sender},
        Arc, Mutex,
    },
    thread,
    time::Duration,
};

use mlsmpm::*;

pub fn run_mlsmpm(tx: Arc<Mutex<Vec<Particle>>>) {
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

    println!("{:?}", settings);

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

    let mut calc = Calculator::new(&settings, Space::new_for_poiseuille(&settings));

    for _i in 0..v_time_steps {
        calc.update();
        *tx.lock().unwrap() = calc.get_particles().to_vec();
    }
}

pub fn run_dambreak(tx: Arc<Mutex<Vec<Particle>>>, step_signal: Option<Receiver<bool>>) {
    let settings = Settings {
        dt: 1e-2,
        gravity: -10.,
        dynamic_viscosity: 1e-6,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 50,
        c: 10.,
        eos_power: 2.,
    };

    println!("{:?}", settings);

    let mut calc = Calculator::new(&settings, Space::new_for_dambreak(&settings));

    loop {
        if let Some(ref sig) = step_signal {
            sig.recv().unwrap();
        }
        calc.update();
        *tx.lock().unwrap() = calc.get_particles().to_vec();
    }
}
