use std::{error::Error, fs, io::Write, path::Path};

use mlsmpm::*;
use nalgebra::Matrix2xX;
use rand::{seq::SliceRandom, Rng, SeedableRng};

fn main() -> Result<(), Box<dyn Error>> {
    const DYNAMIC_VISCOSITY: f64 = 1e-3;
    const DT: f64 = 5e-5;
    let res_list = [50, 100, 250, 500, 1000, 2000];

    for r in res_list.iter().rev() {
        let cell_width = 10. / *r as f64;
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

    let current_time = chrono::Local::now();
    let folder_name = format!(
        "exp_taylorgreen_reset/{:}",
        current_time.format("%Y%m%d_%Hh%Mm%Ss")
    );
    let folder = Path::new(&folder_name);
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let time = 1e-1;

    let pi = std::f64::consts::PI;
    let half_domain_size = 1.;
    fn true_vel(t: f64, x: f64, y: f64, u: f64, pi: f64, nu: f64) -> Vector2<f64> {
        let exp_term = f64::exp(-2. * pi * pi * t / (u * u / nu));
        Vector2::new(
            u * exp_term * f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
            -u * exp_term * f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
        )
    }

    [
        (P2GSchemeType::MLSMPM, G2PSchemeType::MLSMPM),
        (P2GSchemeType::CompactLsmpsLinear, G2PSchemeType::LsmpsFlip1),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::LsmpsFlip1),
        // (P2GSchemeType::Compact3_1, G2PSchemeType::LsmpsFlip1),
        // (P2GSchemeType::Compact1_2, G2PSchemeType::LsmpsFlip2),
        // (P2GSchemeType::Compact2_2, G2PSchemeType::LsmpsFlip2),
        // (P2GSchemeType::Compact3_2, G2PSchemeType::LsmpsFlip2),
        // (
        //     P2GSchemeType::CompactLaplacian2_2,
        //     G2PSchemeType::LsmpsFlip2,
        // ),
        // (
        //     P2GSchemeType::CompactLaplacian3_2,
        //     G2PSchemeType::LsmpsFlip2,
        // ),
        // (P2GSchemeType::CompactV0_2, G2PSchemeType::LsmpsFlip2),
        // (P2GSchemeType::CompactV0_1, G2PSchemeType::LsmpsLinear),
        // (P2GSchemeType::CompactV0_2, G2PSchemeType::LSMPS),
        // (
        //     P2GSchemeType::CompactLsmpsLinear,
        //     G2PSchemeType::LsmpsLinear,
        // ),
        // (P2GSchemeType::Compact1_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::LsmpsLinear),
        // (P2GSchemeType::Compact2_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact3_1, G2PSchemeType::LsmpsLinear),
        // (P2GSchemeType::Compact3_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::CompactLaplacian2_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::CompactLaplacian3_2, G2PSchemeType::LSMPS),
        // // // (P2GSchemeType::MLSMPM, G2PSchemeType::LsmpsLinear),
        // // // // (P2GSchemeType::LsmpsLinear, G2PSchemeType::LsmpsLinear),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::LSMPS),
        // // (P2GSchemeType::CompactLsmpsLinear, G2PSchemeType::LSMPS),
        // // (P2GSchemeType::CompactLsmpsLinear, G2PSchemeType::Lsmps3rd),
        // // // (P2GSchemeType::LSMPS, G2PSchemeType::LSMPS),
        // // // (P2GSchemeType::Lsmps3rd, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Lsmps4th, G2PSchemeType::Lsmps4th),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::CompactLsmps),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::Lsmps3rd),
        // // // (P2GSchemeType::CompactLsmps, G2PSchemeType::MLSMPM),
        // (P2GSchemeType::Compact_0_1, G2PSchemeType::LsmpsLinear),
        // (P2GSchemeType::Compact_0_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_2_2, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_3_2, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_3_3, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_4_3, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_4_4, G2PSchemeType::Lsmps4th),
        // (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::MLSMPM),
        // (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::LSMPS),
        // (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_v_0_1, G2PSchemeType::MLSMPM),
        // (P2GSchemeType::Compact_v_0_1, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_v_3_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_v_3_2, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_v_3_3, G2PSchemeType::Lsmps3rd),
    ]
    .iter()
    .for_each(|&(p2g_scheme, g2p_scheme)| {
        let folder_name = folder_name.clone();
        // thread::spawn(move || {
        let folder = Path::new(&folder_name);
        let results = res_list
            .iter()
            .map(|&grid_width| {
                let folder_name = folder_name.clone();
                // thread::spawn(move || {
                let folder = Path::new(&folder_name);
                let settings = Settings {
                    dt: DT,
                    gravity: 0.,
                    dynamic_viscosity: DYNAMIC_VISCOSITY,
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
                    reset_particle_position: false,
                    ..Default::default()
                };

                println!("{:?}", settings);
                write_settings(folder, p2g_scheme, g2p_scheme, &settings).unwrap();

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

                let linf_error = particles
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
                        .norm()
                    })
                    .fold(0.0 / 0.0, f64::max)
                    / particles
                        .iter()
                        .map(|p| (p.x().x, p.x().y))
                        .map(|(x, y)| {
                            true_vel(time, x, y, half_domain_size, pi, settings.dynamic_viscosity)
                                .norm()
                        })
                        .fold(0.0 / 0.0, f64::max);

                write_final_result(
                    folder,
                    time,
                    half_domain_size,
                    pi,
                    DYNAMIC_VISCOSITY,
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
                                4. <= pos.x && pos.x < 6. && 4. <= pos.y && pos.y < 6.
                            })
                            .collect::<Vec<_>>(),
                    ),
                )
                .unwrap();

                (settings.cell_width() / 2., l2_error, linf_error)
                // })
            })
            .collect::<Vec<_>>();

        write_l2_errors(folder, p2g_scheme, g2p_scheme, results).unwrap();
        // })
    });

    fn write_settings(
        folder: &Path,
        p2g_scheme: P2GSchemeType,
        g2p_scheme: G2PSchemeType,
        settings: &Settings,
    ) -> Result<(), Box<dyn Error>> {
        let mut file = fs::File::create(
            folder.join(format!("settings_{:?}_{:?}.csv", p2g_scheme, g2p_scheme)),
        )?;

        file.write_all(format!("{:?}", settings).as_bytes())?;

        Ok(())
    }

    fn write_l2_errors(
        folder: &Path,
        p2g_scheme: P2GSchemeType,
        g2p_scheme: G2PSchemeType,
        results: Vec<(f64, f64, f64)>,
    ) -> Result<(), Box<dyn Error>> {
        let mut writer = csv::Writer::from_path(
            folder.join(format!("errors_{:?}_{:?}.csv", p2g_scheme, g2p_scheme)),
        )?;
        writer.write_record(&["res", "l2_error", "linf_error"])?;
        for (res, l2_error, linf_error) in results {
            writer.write_record(&[
                res.to_string(),
                l2_error.to_string(),
                linf_error.to_string(),
            ])?;
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

    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(2);

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

            let dvxdxx = -pi * pi / (half_domain_size * half_domain_size)
                * f64::sin(pi * (x - 5.) / half_domain_size)
                * f64::cos(pi * (y - 5.) / half_domain_size);
            let dvxdxy = -pi * pi / (half_domain_size * half_domain_size)
                * f64::cos(pi * (x - 5.) / half_domain_size)
                * f64::sin(pi * (y - 5.) / half_domain_size);
            let dvxdyy = -pi * pi / (half_domain_size * half_domain_size)
                * f64::sin(pi * (x - 5.) / half_domain_size)
                * f64::cos(pi * (y - 5.) / half_domain_size);
            let dvydxx = pi * pi / (half_domain_size * half_domain_size)
                * f64::cos(pi * (x - 5.) / half_domain_size)
                * f64::sin(pi * (y - 5.) / half_domain_size);
            let dvydxy = pi * pi / (half_domain_size * half_domain_size)
                * f64::sin(pi * (x - 5.) / half_domain_size)
                * f64::cos(pi * (y - 5.) / half_domain_size);
            let dvydyy = pi * pi / (half_domain_size * half_domain_size)
                * f64::cos(pi * (x - 5.) / half_domain_size)
                * f64::sin(pi * (y - 5.) / half_domain_size);

            let x_lsmps = Matrix2xX::from_row_slice(&[
                velocity.x, c.m11, c.m12, dvxdxx, dvxdxy, dvxdyy, velocity.y, c.m21, c.m22, dvydxx,
                dvydxy, dvydyy,
            ]);

            let p = Particle::new_with_mass_velocity_c_lsmps(
                Vector2::new(x, y),
                (settings.rho_0 * (half_domain_size * 2.) * (half_domain_size * 2.))
                    / (num_x * num_x) as f64,
                velocity,
                c,
                x_lsmps,
            );
            particles.push(p);
        }
    }
    particles.shuffle(&mut rng);

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
