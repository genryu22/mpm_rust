use std::{thread, vec};

use mlsmpm::*;
use rand::Rng;

fn main() {
    let mut handles = vec![];

    for s in [100, 200, 400, 500, 800] {
        let handle = thread::spawn(move || {
            let settings = Settings {
                dt: 1e-4,
                gravity: 0.,
                dynamic_viscosity: 1e-2,
                alpha: 0.,
                affine: true,
                space_width: 10.,
                grid_width: s,
                rho_0: 1.,
                c: 1e1,
                eos_power: 4.,
                boundary_mirror: false,
                vx_zero: false,
                weight_type: WeightType::QuadraticBSpline,
                p2g_scheme: P2GSchemeType::MLSMPM,
                g2p_scheme: G2PSchemeType::MLSMPM,
            };

            println!("{:?}", settings);

            let v_time_steps = (0.1 / settings.dt) as u32;

            let mut calc = Calculator::new(&settings, new_for_taylor_green(&settings));
            calc.start(v_time_steps);

            let PI = std::f64::consts::PI;
            let half_domain_size = 1.;
            let nu = settings.dynamic_viscosity / settings.rho_0;
            fn true_vel(t: f64, x: f64, y: f64, U: f64, PI: f64, nu: f64) -> Vector2<f64> {
                let exp_term = f64::exp(-2. * PI * PI * t / (U * U / nu));
                Vector2::new(
                    U * exp_term * f64::sin(PI * (x - 5.) / U) * f64::cos(PI * (y - 5.) / U),
                    -U * exp_term * f64::cos(PI * (x - 5.) / U) * f64::sin(PI * (y - 5.) / U),
                )
            }

            file::write_particles(
                calc.get_particles(),
                v_time_steps as usize,
                "exp_taylorgreen",
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

pub fn new_for_taylor_green(settings: &Settings) -> Space {
    let grid_width = settings.grid_width;

    let PI = std::f64::consts::PI;
    let half_domain_size = 1.;

    let pos_x_min = 5. - half_domain_size;
    let pos_x_max = 5. + half_domain_size;
    let num_x = (half_domain_size * 2. / (settings.cell_width() / 2.)) as usize;
    let p_dist = half_domain_size * 2. / (num_x as f64);

    let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

    let mut rng = rand::thread_rng();

    for i_y in 0..num_x {
        for i_x in 0..num_x {
            let mut x = p_dist * (i_x as f64 + 0.5) as f64 + pos_x_min;
            let mut y = p_dist * (i_y as f64 + 0.5) as f64 + pos_x_min;
            if false {
                x += rng.gen_range(-1.0..=1.0) * p_dist * 0.2;
                y += rng.gen_range(-1.0..=1.0) * p_dist * 0.2;
            }
            let velocity = Vector2::new(
                f64::sin(PI * (x - 5.) / half_domain_size)
                    * f64::cos(PI * (y - 5.) / half_domain_size),
                -f64::cos(PI * (x - 5.) / half_domain_size)
                    * f64::sin(PI * (y - 5.) / half_domain_size),
            );

            let c = {
                let k = PI / half_domain_size;
                let c11 = k
                    * f64::cos(PI * (x - 5.) / half_domain_size)
                    * f64::cos(PI * (y - 5.) / half_domain_size);
                let c12 = -k
                    * f64::sin(PI * (x - 5.) / half_domain_size)
                    * f64::sin(PI * (y - 5.) / half_domain_size);
                let c21 = k
                    * f64::sin(PI * (x - 5.) / half_domain_size)
                    * f64::sin(PI * (y - 5.) / half_domain_size);
                let c22 = -k
                    * f64::cos(PI * (x - 5.) / half_domain_size)
                    * f64::cos(PI * (y - 5.) / half_domain_size);

                Matrix2::new(c11, c12, c21, c22)
            };

            let p = Particle::new_with_mass_velocity_c(
                Vector2::new(x, y),
                (settings.rho_0 * (half_domain_size * 2.) * (half_domain_size * 2.))
                    / (num_x * num_x) as f64,
                velocity,
                c,
            );
            particles.push(p);
        }
    }

    let mut grid: Vec<Node> = Vec::with_capacity((grid_width + 1) * (grid_width + 1));
    for i in 0..(grid_width + 1) * (grid_width + 1) {
        grid.push(Node::new((i % (grid_width + 1), i / (grid_width + 1))));
    }

    Space::new(
        grid,
        particles,
        vec![],
        vec![],
        Some(PeriodicBoundaryRect::new(
            pos_x_min, pos_x_max, pos_x_min, pos_x_max,
        )),
        0,
    )
}
