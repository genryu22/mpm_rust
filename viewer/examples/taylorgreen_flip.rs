use mlsmpm::*;

fn main() {
    const DYNAMIC_VISCOSITY: f64 = 1e-2;
    let settings = Settings {
        dt: 1e-4,
        gravity: 0.,
        dynamic_viscosity: DYNAMIC_VISCOSITY,
        alpha: 1.,
        affine: true,
        space_width: 10.,
        grid_width: 100,
        rho_0: 1.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: false,
        vx_zero: false,
        weight_type: WeightType::QuadraticBSpline,
        effect_radius: 2,
        p2g_scheme: P2GSchemeType::MLSMPM,
        g2p_scheme: G2PSchemeType::MLSMPM,
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
        ..Default::default()
    };

    let space = new_for_taylor_green(&settings);
    viewer::run_window_bevy(settings.space_width, settings, space, 500.);
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
