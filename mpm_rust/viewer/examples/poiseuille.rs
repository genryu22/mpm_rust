use mlsmpm::*;

fn main() {
    let settings = Settings {
        dt: 5e-3,
        gravity: -1e-2,
        dynamic_viscosity: 1e-2,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 200,
        rho_0: 1.,
        c: 0.,
        eos_power: 0.,
        boundary_mirror: false,
        vx_zero: true,
        weight_type: WeightType::QuadraticBSpline,
        p2g_scheme: P2GSchemeType::MLSMPM,
        g2p_scheme: G2PSchemeType::MLSMPM,
        pressure: Some(|_, _| 0.),
    };

    let space = new_for_poiseuille(&settings);
    viewer::run_window_bevy(settings.space_width, settings, space, 1000.);
}

pub fn new_for_poiseuille(settings: &Settings) -> Space {
    let grid_width = settings.grid_width;
    let cell_size = settings.cell_width();

    let p_dist = cell_size / 2.;

    let pos_x_min = 4.5;
    let pos_x_max = 5.5;
    let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

    let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

    for i_y in 0..num_x {
        for i_x in 0..num_x {
            let mut p = Particle::new_with_mass(
                Vector2::new(
                    (pos_x_max - pos_x_min) * (i_x as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_x_max - pos_x_min) * (i_y as f64 + 0.5) / num_x as f64 + pos_x_min,
                ),
                (settings.rho_0 * (pos_x_max - pos_x_min) * (pos_x_max - pos_x_min))
                    / (num_x * num_x) as f64,
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
        vec![
            SlipBoundary::new(4.5, Direction::X, true, true, true),
            SlipBoundary::new(5.5, Direction::X, false, true, true),
        ],
        vec![PeriodicBoundary::new(
            BoundaryLine::new(4.5, true),
            BoundaryLine::new(5.5, false),
            Direction::Y,
        )],
        None,
        0,
    )
}
