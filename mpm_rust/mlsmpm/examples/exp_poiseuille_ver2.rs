use std::{error::Error, fs, path::Path};

use mlsmpm::*;
use rand::Rng;
use rayon::prelude::*;

fn main() -> Result<(), Box<dyn Error>> {
    let folder = Path::new("exp_poiseuille_ver2");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let g = -1e-2;
    let dynamic_viscosity = 1e-2;

    fn true_vel(x: f64, y: f64, g: f64, nu: f64) -> Vector2<f64> {
        let max = 1.5 * g * 1. * 1. / ((nu / 1.) * 12.);

        Vector2::new(0., -max / 0.25 * (x - 5.).powi(2) + max)
    }

    let result = [
        (P2GSchemeType::MLSMPM, G2PSchemeType::MLSMPM),
        // (P2GSchemeType::MLSMPM, G2PSchemeType::LsmpsLinear),
        // (P2GSchemeType::LsmpsLinear, G2PSchemeType::LsmpsLinear),
        // (
        //     P2GSchemeType::CompactLsmpsLinear,
        //     G2PSchemeType::LsmpsLinear,
        // ),
        // (
        //     P2GSchemeType::CompactOnlyVelocity,
        //     G2PSchemeType::LsmpsLinear,
        // ),
        (P2GSchemeType::LSMPS, G2PSchemeType::LSMPS),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::LSMPS),
    ]
    .par_iter()
    .map(|&(p2g_scheme, g2p_scheme)| {
        let result = [100, 200, 400, 500, 800]
            .par_iter()
            .map(|&grid_width| {
                let settings = Settings {
                    dt: 5e-3,
                    gravity: g,
                    dynamic_viscosity,
                    alpha: 0.,
                    affine: true,
                    space_width: 10.,
                    grid_width,
                    rho_0: 1.,
                    c: 0.,
                    eos_power: 0.,
                    boundary_mirror: true,
                    vx_zero: true,
                    weight_type: WeightType::QuadraticBSpline,
                    effect_radius: 3,
                    p2g_scheme,
                    g2p_scheme,
                    pressure: Some(|_, _| 0.),
                    ..Default::default()
                };

                println!("{:?}", settings);

                let v_time_steps = (1. / settings.dynamic_viscosity / settings.dt).ceil() as u32;

                // println!("粘性時間: L^2/mu = {} steps", v_time_steps);
                // println!(
                //     "平均流速の理論解: U_mean = g*L*L/(nu*12) = {}",
                //     settings.gravity * 1. * 1. / ((settings.dynamic_viscosity / 1.) * 12.)
                // );
                // println!(
                //     "最大流速の理論解: U_max = 1.5*U_mean = {}",
                //     1.5 * settings.gravity * 1. * 1. / ((settings.dynamic_viscosity / 1.) * 12.)
                // );

                let space = new_for_poiseuille(&settings);
                let mut calc = Calculator::new(&settings, space);
                calc.start(v_time_steps);

                let particles = calc.get_particles();
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

                (
                    settings.cell_width() / 2.,
                    l2_error,
                    calc.get_grid()
                        .iter()
                        .enumerate()
                        .map(|(index, node)| {
                            (
                                Vector2::<f64>::new(
                                    (index % (settings.grid_width + 1)) as f64,
                                    (index / (settings.grid_width + 1)) as f64,
                                ) * settings.cell_width(),
                                node.clone(),
                            )
                        })
                        .filter(|(pos, _)| {
                            4.5 <= pos.x && pos.x <= 5.5 && 4.5 <= pos.y && pos.y <= 5.5
                        })
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>();

        (p2g_scheme, g2p_scheme, result)
    })
    .collect::<Vec<_>>();

    for (p2g_scheme, g2p_shceme, result) in result {
        let mut writer =
            csv::Writer::from_path(folder.join(format!("{:?}_{:?}.csv", p2g_scheme, g2p_shceme)))?;
        writer.write_record(&["res", "l2_error"])?;
        for (index, (res, l2_error, grid)) in result.iter().enumerate() {
            writer.write_record(&[res.to_string(), l2_error.to_string()])?;

            if index == result.len() - 1 {
                let mut writer = csv::Writer::from_path(folder.join(format!(
                    "final_result_{:?}_{:?}.csv",
                    p2g_scheme, g2p_shceme
                )))?;
                writer.write_record(&["x", "y", "vx", "vy", "t_vx", "t_vy"])?;
                for (pos, n) in grid {
                    let true_vel = true_vel(pos.x, pos.y, g, dynamic_viscosity);
                    writer.write_record(&[
                        pos.x.to_string(),
                        pos.y.to_string(),
                        n.v.x.to_string(),
                        n.v.y.to_string(),
                        true_vel.x.to_string(),
                        true_vel.y.to_string(),
                    ])?;
                }
                writer.flush()?;
            }
        }
        writer.flush()?;
    }

    Ok(())
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
