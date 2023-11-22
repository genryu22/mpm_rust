use std::{error::Error, fs, path::Path, vec};

use mlsmpm::*;

fn main() -> Result<(), Box<dyn Error>> {
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
        weight_type: WeightType::QuadraticBSpline,
        effect_radius: 3,
        p2g_scheme: P2GSchemeType::MLSMPM,
        g2p_scheme: G2PSchemeType::MLSMPM,
        pressure: None,
        ..Default::default()
    };

    println!("{:?}", settings);

    let folder = Path::new("exp_dambreak");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let mut particle_writer = csv::Writer::from_path(folder.join("dambreak.csv"))?;
    particle_writer.write_record(&["nt", "z"])?;

    let mut calc = Calculator::new(&settings, new_for_dambreak_experiment(&settings));

    particle_writer.write_record(&["0", &(calc.get_max_x() - 3.).to_string()])?;

    for i in 0..7900 {
        calc.update();

        let step = (i + 1) as f64;
        particle_writer.write_record(&[
            &(step * settings.dt * f64::sqrt(-2. * settings.gravity)).to_string(),
            &(calc.get_max_x() - 3.).to_string(),
        ])?;
    }

    particle_writer.flush()?;

    Ok(())
}

pub fn new_for_dambreak_experiment(settings: &Settings) -> Space {
    let grid_width = settings.grid_width;
    let cell_size = settings.cell_width();

    let p_dist = cell_size / 2.;

    let pos_x_min = 3.0;
    let pos_x_max = 4.0;
    let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

    let pos_y_min = 4.0;
    let pos_y_max = 6.0;
    let num_y = ((pos_y_max - pos_y_min) / p_dist + 0.5) as usize;

    let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

    for i_y in 0..num_y {
        for i_x in 0..num_x {
            let p = Particle::new_with_mass(
                Vector2::new(
                    (pos_x_max - pos_x_min) * (i_x as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_y_max - pos_y_min) * (i_y as f64 + 0.5) / num_y as f64 + pos_y_min,
                ),
                (settings.rho_0 * (pos_x_max - pos_x_min) * (pos_y_max - pos_y_min))
                    / (num_x * num_y) as f64,
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
            SlipBoundary::new(3., Direction::X, true, true, false),
            SlipBoundary::new(7., Direction::X, false, true, false),
            SlipBoundary::new(4., Direction::Y, true, false, false),
            SlipBoundary::new(6.5, Direction::Y, false, false, false),
        ],
        vec![],
        None,
        0,
    )
}
