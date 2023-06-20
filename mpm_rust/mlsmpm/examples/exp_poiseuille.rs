use std::{thread, vec};

use mlsmpm::*;

fn main() {
    let mut handles = vec![];

    for s in [100, 200, 400, 500, 800] {
        let handle = thread::spawn(move || {
            let settings = Settings {
                dt: 5e-3,
                gravity: -1e-2,
                dynamic_viscosity: 1e-2,
                alpha: 0.,
                affine: true,
                space_width: 10.,
                grid_width: s as usize,
                rho_0: 1.,
                c: 0.,
                eos_power: 0.,
                boundary_mirror: true,
                vx_zero: true,
                weight_type: WeightType::QuadraticBSpline,
                scheme: SchemeType::MLSMPM,
            };

            println!("{:?}", settings);

            let v_time_steps = (1. / settings.dynamic_viscosity / settings.dt).ceil() as u32;
            println!("粘性時間: L^2/mu = {} steps", v_time_steps);

            let mut calc = Calculator::new(&settings, Space::new_for_poiseuille(&settings));
            calc.start(v_time_steps);

            file::write_particles(
                calc.get_particles(),
                v_time_steps as usize,
                "exp_poiseuille",
                &format!("{}", s),
            )
            .unwrap();
        });

        handles.push(handle);
    }

    for h in handles.into_iter() {
        h.join().unwrap();
    }
}
