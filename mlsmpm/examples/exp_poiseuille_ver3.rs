use std::{
    error::Error,
    fs,
    path::{Path, PathBuf},
};

use mlsmpm::*;

fn main() -> Result<(), Box<dyn Error>> {
    let g = -1e-2;
    let dynamic_viscosity = 1e-2;

    let schemes = [
        (P2GSchemeType::MLSMPM, G2PSchemeType::MLSMPM),
        (P2GSchemeType::MLSMPM, G2PSchemeType::LsmpsLinear),
        (P2GSchemeType::LsmpsLinear, G2PSchemeType::LsmpsLinear),
        (P2GSchemeType::LSMPS, G2PSchemeType::LSMPS),
        (P2GSchemeType::Lsmps3rd, G2PSchemeType::LSMPS),
    ];

    let res_list = [100, 200, 400, 500, 800];
    let dt_list = res_list
        .iter()
        .map(|r| {
            let u_max = -g * 1. * 1. / 8. / dynamic_viscosity;
            let cell_width = 10. / *r as f64;
            let dx = cell_width / 2.;

            f64::min(dx / 2. / u_max, dx * dx / 10. / dynamic_viscosity)
        })
        .collect::<Vec<_>>();

    // for r in res_list.iter().rev() {
    //     let u_max = -g * 1. * 1. / 8. / dynamic_viscosity;
    //     let cell_width = 10. / *r as f64;
    //     let dx = cell_width / 2.;
    //     let dt_max = f64::min(dx / 2. / u_max, dx * dx / 10. / dynamic_viscosity);
    //     assert!(dt <= dt_max, "dt = {} > {}", dt, dt_max);
    // }

    let folder = Path::new("exp_poiseuille_ver3");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }
    let current_time = chrono::Local::now();
    let folder = folder.join(current_time.format("%Y%m%d_%Hh%Mm%Ss").to_string());
    if !folder.exists() {
        fs::create_dir(folder.clone())?;
    }

    let result = res_list
        .iter()
        .enumerate()
        .map(|(i, &grid_width)| {
            let dt = dt_list[i];
            schemes.clone().map(|(p2g_scheme, g2p_scheme)| {
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

                let v_time_steps = (1. / settings.dynamic_viscosity / settings.dt).ceil() as u32;

                let space = new_for_poiseuille(&settings);
                let mut calc = Calculator::new(&settings, space);
                calc.start(v_time_steps);

                let particles = calc.get_particles();

                write_final_result(
                    &settings, &folder, p2g_scheme, g2p_scheme, grid_width, particles,
                )
                .unwrap();

                let l2_error = f64::sqrt(
                    particles
                        .iter()
                        .map(|p| {
                            let x = p.x().x;
                            let y = p.x().y;
                            (p.v() - true_vel(x, y, settings.gravity, settings.dynamic_viscosity))
                                .norm_squared()
                        })
                        .sum::<f64>()
                        / particles
                            .iter()
                            .map(|p| (p.x().x, p.x().y))
                            .map(|(x, y)| {
                                true_vel(x, y, settings.gravity, settings.dynamic_viscosity)
                                    .norm_squared()
                            })
                            .sum::<f64>(),
                );

                let linf_error = f64::sqrt(
                    particles
                        .iter()
                        .map(|p| {
                            let x = p.x().x;
                            let y = p.x().y;
                            (p.v() - true_vel(x, y, settings.gravity, settings.dynamic_viscosity))
                                .norm()
                        })
                        .fold(0.0 / 0.0, f64::max)
                        / particles
                            .iter()
                            .map(|p| (p.x().x, p.x().y))
                            .map(|(x, y)| {
                                true_vel(x, y, settings.gravity, settings.dynamic_viscosity).norm()
                            })
                            .fold(0.0 / 0.0, f64::max),
                );

                (l2_error, linf_error)
            })
        })
        .collect::<Vec<_>>();

    {
        let file = folder.join("l2.csv");
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(
            &([
                vec!["res".to_string()],
                schemes
                    .map(|(p2g, g2p)| format!("{:?}_{:?}", p2g, g2p))
                    .to_vec(),
            ]
            .concat()),
        )?;
        for i in 0..res_list.len() {
            writer.write_record(
                &([
                    vec![res_list[i].to_string()],
                    result[i].map(|(l2_error, _)| l2_error.to_string()).to_vec(),
                ]
                .concat()),
            )?;
        }
        writer.flush()?;
    }

    {
        let file = folder.join("linf.csv");
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(
            &([
                vec!["res".to_string()],
                schemes
                    .map(|(p2g, g2p)| format!("{:?}_{:?}", p2g, g2p))
                    .to_vec(),
            ]
            .concat()),
        )?;
        for i in 0..res_list.len() {
            writer.write_record(
                &([
                    vec![res_list[i].to_string()],
                    result[i]
                        .map(|(_, linf_error)| linf_error.to_string())
                        .to_vec(),
                ]
                .concat()),
            )?;
        }
        writer.flush()?;
    }

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
