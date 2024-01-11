use std::vec;

use mlsmpm::*;
use rand::Rng;

fn main() {
    const DYNAMIC_VISCOSITY: f64 = 1e-3;
    const DT: f64 = 5e-4;
    const GRID_WIDTH: usize = 1000;

    {
        let cell_width = 10. / GRID_WIDTH as f64;
        assert!(
            DT <= f64::min(
                (cell_width / 2.) / 2. / 1.,
                (cell_width / 2.).powi(2) / 10. / DYNAMIC_VISCOSITY
            ),
            "dt = {} > {}",
            DT,
            f64::min(
                (cell_width / 2.) / 2. / 1.,
                (cell_width / 2.).powi(2) / 10. / DYNAMIC_VISCOSITY
            )
        );
    }

    let settings = Settings {
        dt: DT,
        gravity: 0.,
        dynamic_viscosity: DYNAMIC_VISCOSITY,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: GRID_WIDTH,
        rho_0: 1.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: false,
        vx_zero: false,
        weight_type: WeightType::QuinticBSpline,
        effect_radius: 5,
        p2g_scheme: P2GSchemeType::LSMPS,
        g2p_scheme: G2PSchemeType::LSMPS,
        pressure: Some(|p, time| {
            let pi = std::f64::consts::PI;
            let l = 1.;
            let rho = 1.;
            let u = 1.;

            let (x, y) = (p.x().x - 5., p.x().y - 5.);

            rho * u * u / 4.
                * f64::exp(-4. * pi * pi * time * DYNAMIC_VISCOSITY / (l * l))
                * (f64::cos(2. * pi * x / l) + f64::cos(2. * pi * y / l))
        }),
        pressure_grad: Some(|x, y, time| {
            let pi = std::f64::consts::PI;
            let l = 1.;
            let rho = 1.;
            let u = 1.;

            let (x, y) = (x - 5., y - 5.);

            let p_dx = rho * u * u * pi / 2. / l
                * f64::exp(-4. * pi * pi * time * DYNAMIC_VISCOSITY / (l * l))
                * (-f64::sin(2. * pi * x / l));
            let p_dy = rho * u * u * pi / 2. / l
                * f64::exp(-4. * pi * pi * time * DYNAMIC_VISCOSITY / (l * l))
                * (-f64::sin(2. * pi * y / l));

            Vector2::new(p_dx, p_dy)
        }),
        parallel: false,
        ..Default::default()
    };

    println!("{:?}", settings);

    let v_time_steps = (0.1 / settings.dt) as u32;

    let mut calc = Calculator::new(&settings, new_for_taylor_green(&settings));
    calc.start(v_time_steps);
}

pub fn new_for_taylor_green(settings: &Settings) -> Space {
    let grid_width = settings.grid_width;

    let pi = std::f64::consts::PI;
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
            if true {
                x += rng.gen_range(-1.0..=1.0) * p_dist * 0.2;
                y += rng.gen_range(-1.0..=1.0) * p_dist * 0.2;
            }
            let velocity = Vector2::new(
                f64::sin(pi * (x - 5.) / half_domain_size)
                    * f64::cos(pi * (y - 5.) / half_domain_size),
                -f64::cos(pi * (x - 5.) / half_domain_size)
                    * f64::sin(pi * (y - 5.) / half_domain_size),
            );

            let c = {
                let k = pi / half_domain_size;
                let c11 = k
                    * f64::cos(pi * (x - 5.) / half_domain_size)
                    * f64::cos(pi * (y - 5.) / half_domain_size);
                let c12 = -k
                    * f64::sin(pi * (x - 5.) / half_domain_size)
                    * f64::sin(pi * (y - 5.) / half_domain_size);
                let c21 = k
                    * f64::sin(pi * (x - 5.) / half_domain_size)
                    * f64::sin(pi * (y - 5.) / half_domain_size);
                let c22 = -k
                    * f64::cos(pi * (x - 5.) / half_domain_size)
                    * f64::cos(pi * (y - 5.) / half_domain_size);

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
