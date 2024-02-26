use std::{
    error::Error,
    fs,
    path::{Path, PathBuf},
};

use mlsmpm::*;

fn main() -> Result<(), Box<dyn Error>> {
    let g = -1e-2;
    let dynamic_viscosity = 1e-2;

    let schemes = (P2GSchemeType::MLSMPM, G2PSchemeType::MLSMPM);

    let grid_width = 500;
    let dt = {
        let u_max = -g * 1. * 1. / 8. / dynamic_viscosity;
        let cell_width = 10. / grid_width as f64;
        let dx = cell_width / 2.;

        f64::min(dx / 2. / u_max, dx * dx / 10. / dynamic_viscosity)
    };

    let folder = Path::new("exp_poiseuille_one");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }
    let current_time = chrono::Local::now();
    let folder = folder.join(current_time.format("%Y%m%d_%Hh%Mm%Ss").to_string());
    if !folder.exists() {
        fs::create_dir(folder.clone())?;
    }

    let (p2g_scheme, g2p_scheme) = schemes;
    let settings = Settings {
        dt,
        gravity: g,
        dynamic_viscosity,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width,
        rho_0: 1.,
        c: 0.,
        eos_power: 0.,
        boundary_mirror: false,
        vx_zero: true,
        weight_type: WeightType::QuadraticBSpline,
        effect_radius: 2,
        p2g_scheme,
        g2p_scheme,
        pressure: Some(|_, _| 0.),
        ..Default::default()
    };

    println!("{:?}_{:?}_{}_{}s", p2g_scheme, g2p_scheme, grid_width, dt);

    let v_time_steps = 1000; //(1. / settings.dynamic_viscosity / settings.dt).ceil() as u32;

    let space = new_for_poiseuille(&settings);
    let mut calc = Calculator::new(&settings, space);
    calc.start(v_time_steps);

    let particles = calc.get_particles();

    write_final_result(
        &settings, &folder, p2g_scheme, g2p_scheme, grid_width, particles,
    )?;

    fn true_vel(x: f64, _y: f64, g: f64, nu: f64) -> Vector2<f64> {
        let max = 1.5 * g * 1. * 1. / ((nu / 1.) * 12.);

        Vector2::new(0., -max / 0.25 * (x - 5.).powi(2) + max)
    }

    fn write_final_result(
        settings: &Settings,
        folder: &PathBuf,
        p2g: P2GSchemeType,
        g2p: G2PSchemeType,
        res: usize,
        particles: &Vec<Particle>,
    ) -> Result<(), Box<dyn Error>> {
        let folder = folder.join("distribution");
        if !folder.exists() {
            fs::create_dir(folder.clone())?;
        }
        let folder = folder.join(format!("{:?}_{:?}", p2g, g2p));
        if !folder.exists() {
            fs::create_dir(folder.clone())?;
        }
        let file = folder.join(format!("{}.csv", res));
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(&["x", "y", "u", "v", "u_exact", "v_exact"])?;
        for p in particles.iter() {
            let exact = true_vel(
                p.x().x,
                p.x().y,
                settings.gravity,
                settings.dynamic_viscosity,
            );
            writer.write_record(&[
                p.x().x.to_string(),
                p.x().y.to_string(),
                p.v().x.to_string(),
                p.v().y.to_string(),
                exact.x.to_string(),
                exact.y.to_string(),
            ])?;
        }
        writer.flush()?;

        Ok(())
    }

    Ok(())
}

fn new_for_poiseuille(settings: &Settings) -> Space {
    let grid_width = settings.grid_width;
    let cell_size = settings.cell_width();

    let p_dist = cell_size / 2.;

    let pos_x_min = 4.5;
    let pos_x_max = 5.5;
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
