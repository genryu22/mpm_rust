use std::{error::Error, fs, path::Path};

use mlsmpm::*;
use rand::Rng;
use rayon::prelude::*;

fn main() -> Result<(), Box<dyn Error>> {
    let folder = Path::new("exp_taylorgreen_ver2");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let result = [
        (P2GSchemeType::MLSMPM, G2PSchemeType::MLSMPM),
        (P2GSchemeType::MLSMPM, G2PSchemeType::LsmpsLinear),
        (P2GSchemeType::LsmpsLinear, G2PSchemeType::LsmpsLinear),
        (
            P2GSchemeType::CompactLsmpsLinear,
            G2PSchemeType::LsmpsLinear,
        ),
        (P2GSchemeType::LSMPS, G2PSchemeType::LSMPS),
        (P2GSchemeType::CompactLsmps, G2PSchemeType::LSMPS),
        (P2GSchemeType::CompactLsmps, G2PSchemeType::CompactLsmps),
    ]
    .par_iter()
    .map(|&(p2g_scheme, g2p_scheme)| {
        let result = [25, 50, 100] //, 200, 250, 400, 500, 800]
            .par_iter()
            .map(|&grid_width| {
                let settings = Settings {
                    dt: 1e-4,
                    gravity: 0.,
                    dynamic_viscosity: 1e-2,
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
                    p2g_scheme,
                    g2p_scheme,
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

                println!("{:?}", settings);

                let time = 1e-1;
                let v_time_steps = (time / settings.dt) as u32;

                let space = new_for_taylor_green(&settings);
                let mut calc = Calculator::new(&settings, space);
                calc.start(v_time_steps);

                let PI = std::f64::consts::PI;
                let half_domain_size = 1.;
                fn true_vel(t: f64, x: f64, y: f64, U: f64, PI: f64, nu: f64) -> Vector2<f64> {
                    let exp_term = f64::exp(-2. * PI * PI * t / (U * U / nu));
                    Vector2::new(
                        U * exp_term * f64::sin(PI * (x - 5.) / U) * f64::cos(PI * (y - 5.) / U),
                        -U * exp_term * f64::cos(PI * (x - 5.) / U) * f64::sin(PI * (y - 5.) / U),
                    )
                }

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
                                    PI,
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
                                    PI,
                                    settings.dynamic_viscosity,
                                )
                                .norm_squared()
                            })
                            .sum::<f64>(),
                );

                (settings.cell_width() / 2., l2_error)
            })
            .collect::<Vec<_>>();

        (p2g_scheme, g2p_scheme, result)
    })
    .collect::<Vec<_>>();

    for (p2g_scheme, g2p_shceme, result) in result {
        let mut writer =
            csv::Writer::from_path(folder.join(format!("{:?}_{:?}.csv", p2g_scheme, g2p_shceme)))?;
        writer.write_record(&["res", "l2_error"])?;
        for (res, l2_error) in result {
            writer.write_record(&[res.to_string(), l2_error.to_string()])?;
        }
        writer.flush()?;
    }

    Ok(())
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
                f64::sin(PI * (x - 5.) / half_domain_size)
                    * f64::cos(PI * (y - 5.) / half_domain_size),
                -f64::cos(PI * (x - 5.) / half_domain_size)
                    * f64::sin(PI * (y - 5.) / half_domain_size),
            );

            let c = {
                let k = PI / half_domain_size;
                let c11 = k
                    * f64::cos(PI * (x - 5.) / half_domain_size)
                    * f64::cos(PI * (y - 5.) / half_domain_size);
                let c12 = -k
                    * f64::sin(PI * (x - 5.) / half_domain_size)
                    * f64::sin(PI * (y - 5.) / half_domain_size);
                let c21 = k
                    * f64::sin(PI * (x - 5.) / half_domain_size)
                    * f64::sin(PI * (y - 5.) / half_domain_size);
                let c22 = -k
                    * f64::cos(PI * (x - 5.) / half_domain_size)
                    * f64::cos(PI * (y - 5.) / half_domain_size);

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
