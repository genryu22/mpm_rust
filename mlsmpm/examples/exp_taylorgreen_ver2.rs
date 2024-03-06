use std::{error::Error, fs, path::Path};

use mlsmpm::*;
use rand::Rng;

fn main() -> Result<(), Box<dyn Error>> {
    let folder = Path::new("exp_taylorgreen_ver2");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let time = 1e-2;

    let pi = std::f64::consts::PI;
    let half_domain_size = 1.;
    let dynamic_viscosity = 1e-2;
    fn true_vel(t: f64, x: f64, y: f64, u: f64, pi: f64, nu: f64) -> Vector2<f64> {
        let exp_term = f64::exp(-2. * pi * pi * t / (u * u / nu));
        Vector2::new(
            u * exp_term * f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
            -u * exp_term * f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
        )
    }

    [
        (P2GSchemeType::MLSMPM, G2PSchemeType::MLSMPM),
        (P2GSchemeType::MLSMPM, G2PSchemeType::LsmpsLinear),
        (P2GSchemeType::LsmpsLinear, G2PSchemeType::LsmpsLinear),
        (
            P2GSchemeType::CompactLsmpsLinear,
            G2PSchemeType::LsmpsLinear,
        ),
        (P2GSchemeType::CompactLsmpsLinear, G2PSchemeType::LSMPS),
        (P2GSchemeType::LSMPS, G2PSchemeType::LSMPS),
        (P2GSchemeType::Lsmps3rd, G2PSchemeType::Lsmps3rd),
        (P2GSchemeType::Lsmps4th, G2PSchemeType::Lsmps4th),
        (P2GSchemeType::CompactLsmps, G2PSchemeType::LSMPS),
        (P2GSchemeType::CompactLsmps, G2PSchemeType::CompactLsmps),
        (P2GSchemeType::CompactLsmps, G2PSchemeType::Lsmps3rd),
        (P2GSchemeType::CompactLsmps, G2PSchemeType::MLSMPM),
        (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::MLSMPM),
        (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::LSMPS),
        (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::Lsmps3rd),
    ]
    .iter()
    .for_each(|&(p2g_scheme, g2p_scheme)| {
        let results = [250, 500, 1000]
            .iter()
            .map(|&grid_width| {
                let settings = Settings {
                    dt: 1e-4,
                    gravity: 0.,
                    dynamic_viscosity,
                    alpha: 0.,
                    affine: true,
                    space_width: 10.,
                    grid_width,
                    rho_0: 1.,
                    c: 1e1,
                    eos_power: 4.,
                    boundary_mirror: false,
                    vx_zero: false,
                    weight_type: WeightType::QuadraticBSpline,
                    effect_radius: 2,
                    p2g_scheme,
                    g2p_scheme,
                    pressure: Some(|p, time| {
                        let pi = std::f64::consts::PI;
                        let l = 1.;
                        let rho = 1.;
                        let u = 1.;
                        let nu = 1e-2;

                        let (x, y) = (p.x().x - 5., p.x().y - 5.);

                        rho * u * u / 4.
                            * f64::exp(-4. * pi * pi * time * nu / (l * l))
                            * (f64::cos(2. * pi * x / l) + f64::cos(2. * pi * y / l))
                    }),
                    ..Default::default()
                };

                println!("{:?}", settings);
                let v_time_steps = (time / settings.dt) as u32;

                let space = new_for_taylor_green(&settings);
                let mut calc = Calculator::new(&settings, space);
                calc.start(v_time_steps);

                let particles = calc.get_particles();
                let l2_error = f64::sqrt(
                    particles
                        .iter()
                        .map(|p| {
                            let x = p.x().x;
                            let y = p.x().y;
                            (p.v()
                                - true_vel(
                                    time,
                                    x,
                                    y,
                                    half_domain_size,
                                    pi,
                                    settings.dynamic_viscosity,
                                ))
                            .norm_squared()
                        })
                        .sum::<f64>()
                        / particles
                            .iter()
                            .map(|p| (p.x().x, p.x().y))
                            .map(|(x, y)| {
                                true_vel(
                                    time,
                                    x,
                                    y,
                                    half_domain_size,
                                    pi,
                                    settings.dynamic_viscosity,
                                )
                                .norm_squared()
                            })
                            .sum::<f64>(),
                );

                write_final_result(
                    folder,
                    time,
                    half_domain_size,
                    pi,
                    dynamic_viscosity,
                    p2g_scheme,
                    g2p_scheme,
                    (
                        grid_width,
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
                                4. <= pos.x && pos.x <= 6. && 4. <= pos.y && pos.y <= 6.
                            })
                            .collect::<Vec<_>>(),
                    ),
                )
                .unwrap();

                (settings.cell_width() / 2., l2_error)
            })
            .collect::<Vec<_>>();

        write_l2_errors(folder, p2g_scheme, g2p_scheme, results).unwrap();
    });

    fn write_l2_errors(
        folder: &Path,
        p2g_scheme: P2GSchemeType,
        g2p_scheme: G2PSchemeType,
        results: Vec<(f64, f64)>,
    ) -> Result<(), Box<dyn Error>> {
        let mut writer = csv::Writer::from_path(
            folder.join(format!("errors_{:?}_{:?}.csv", p2g_scheme, g2p_scheme)),
        )?;
        writer.write_record(&["res", "l2_error"])?;
        for (res, l2_error) in results {
            writer.write_record(&[res.to_string(), l2_error.to_string()])?;
        }
        writer.flush()?;

        Ok(())
    }

    fn write_final_result(
        folder: &Path,
        time: f64,
        half_domain_size: f64,
        pi: f64,
        dynamic_viscosity: f64,
        p2g_scheme: P2GSchemeType,
        g2p_shceme: G2PSchemeType,
        result: (usize, Vec<(Vector2<f64>, Node)>),
    ) -> Result<(), Box<dyn Error>> {
        let (res, grid) = result;

        let mut writer = csv::Writer::from_path(folder.join(format!(
            "final_dist_{:?}_{:?}_{:}.csv",
            p2g_scheme, g2p_shceme, res
        )))?;
        writer.write_record(&["x", "y", "vx", "vy", "t_vx", "t_vy"])?;
        for (pos, n) in grid {
            let true_vel = true_vel(time, pos.x, pos.y, half_domain_size, pi, dynamic_viscosity);
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

        Ok(())
    }

    Ok(())
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
            if false {
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