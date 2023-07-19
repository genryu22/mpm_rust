use mlsmpm::*;

fn main() {
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
        weight_type: WeightType::CubicBSpline,
        p2g_scheme: P2GSchemeType::MLSMPM,
        g2p_scheme: G2PSchemeType::MLSMPM,
        pressure: None,
    };

    let space = new_for_dambreak(&settings);
    viewer::run_window_bevy(settings.space_width, settings, space, 100.);
}

pub fn new_for_dambreak(settings: &Settings) -> Space {
    let grid_width = settings.grid_width;
    let cell_size = settings.cell_width();

    let p_dist = cell_size / 2.;

    let pos_x_min = 1.1;
    let pos_x_max = 5.1;
    let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

    let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

    for i_y in 0..num_x {
        for i_x in 0..num_x {
            let p = Particle::new_with_mass(
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
            SlipBoundary::new(1., Direction::X, true, true, false),
            SlipBoundary::new(9., Direction::X, false, true, false),
            SlipBoundary::new(1., Direction::Y, true, true, false),
            SlipBoundary::new(9., Direction::Y, false, true, false),
        ],
        vec![],
        None,
        0,
    )
}
