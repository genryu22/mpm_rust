use mlsmpm::*;

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
        p2g_scheme: P2GSchemeType::LsmpsLinear,
        g2p_scheme: G2PSchemeType::LsmpsLinear,
        pressure: Some(|p, time| {
            let PI = std::f64::consts::PI;
            let L = 1.;
            let rho = 1.;
            let U = 1.;
            let nu = 1e-2;

            let (x, y) = (p.x().x - 5., p.x().y - 5.);

            rho * U * U / 4.
                * f64::exp(-4. * PI * PI * time * nu / (L * L))
                * (f64::cos(2. * PI * x / L) + f64::cos(2. * PI * y / L))
        }),
    };

    let space = new_for_taylor_green(&settings);
    viewer::run_window_bevy(settings.space_width, settings, space, 450.);
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
