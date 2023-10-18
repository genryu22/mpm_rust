use std::{sync::mpsc, thread};

use mlsmpm::*;
use viewer::{
    executor::{Executor, StepExecutor},
    window_bevy, Snapshot,
};

fn main() {
    let settings = Settings {
        dt: 1e-4,
        gravity: 0.,
        dynamic_viscosity: 1e-2,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 100,
        rho_0: 1.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: false,
        vx_zero: false,
        weight_type: WeightType::QuadraticBSpline,
        effect_radius: 3,
        p2g_scheme: P2GSchemeType::MLSMPM,
        g2p_scheme: G2PSchemeType::MLSMPM,
        pressure: None,
        ..Default::default()
    };

    let space = new_for_taylor_green(&settings);

    let count = space.get_particle_count() as u32;

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
        let step_executor = ComparisonExecutor::new(calc, &settings);
        let mut executor = Executor::new(data_sender, step_receiver);
        executor.start(step_executor);
    });

    window_bevy::run(data_receiver, Some(step_sender), |snapshot| {}, 250.);
}

pub struct ComparisonExecutor<'a, 'b> {
    steps: usize,
    calculator: Calculator<'a>,
    settings: &'b Settings,
}

impl<'a, 'b> ComparisonExecutor<'a, 'b> {
    pub fn new(calculator: Calculator<'a>, settings: &'b Settings) -> Self {
        Self {
            steps: 0,
            calculator,
            settings,
        }
    }
}

impl<'a, 'b> StepExecutor for ComparisonExecutor<'a, 'b> {
    fn step(&mut self) -> Snapshot {
        self.calculator.update();
        self.steps += 1;

        let PI = std::f64::consts::PI;
        let half_domain_size = 1.;
        let dynamic_viscosity = 1e-2;
        fn true_vel(t: f64, x: f64, y: f64, U: f64, PI: f64, nu: f64) -> Vector2<f64> {
            let exp_term = f64::exp(-2. * PI * PI * t / (U * U / nu));
            Vector2::new(
                U * exp_term * f64::sin(PI * (x - 5.) / U) * f64::cos(PI * (y - 5.) / U),
                -U * exp_term * f64::cos(PI * (x - 5.) / U) * f64::sin(PI * (y - 5.) / U),
            )
        }

        let particles = self
            .calculator
            .get_particles()
            .iter()
            .map(|p| {
                let pos = p.x().clone();
                Particle::new_with_mass_velocity(
                    pos - Vector2::<f64>::identity() * 1.2,
                    1.,
                    true_vel(
                        self.steps as f64 * self.settings.dt,
                        pos.x,
                        pos.y,
                        half_domain_size,
                        PI,
                        dynamic_viscosity,
                    ),
                )
            })
            .chain(self.calculator.get_particles().iter().map(|p| {
                let pos = p.x().clone();
                Particle::new_with_mass_velocity(
                    pos + Vector2::<f64>::identity() * 1.2,
                    1.,
                    p.v().clone(),
                )
            }))
            .collect::<Vec<_>>();

        Snapshot::new(particles, self.calculator.get_grid().to_vec(), self.steps)
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

    for i_y in 0..num_x {
        for i_x in 0..num_x {
            let x = p_dist * (i_x as f64 + 0.5) as f64 + pos_x_min;
            let y = p_dist * (i_y as f64 + 0.5) as f64 + pos_x_min;
            let p = Particle::new_with_mass_velocity(
                Vector2::new(x, y),
                (settings.rho_0 * (half_domain_size * 2.) * (half_domain_size * 2.))
                    / (num_x * num_x) as f64,
                Vector2::new(
                    f64::sin(PI * (x - 5.) / half_domain_size)
                        * f64::cos(PI * (y - 5.) / half_domain_size),
                    -f64::cos(PI * (x - 5.) / half_domain_size)
                        * f64::sin(PI * (y - 5.) / half_domain_size),
                ),
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
