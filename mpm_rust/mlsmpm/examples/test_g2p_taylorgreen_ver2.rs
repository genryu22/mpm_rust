use std::{error::Error, fs, path::Path};

use mlsmpm::*;
use rand::Rng;

fn main() -> Result<(), Box<dyn Error>> {
    let current_time = chrono::Local::now();
    let folder_name = "exp_g2p_taylorgreen_ver2";
    let folder = Path::new(&folder_name);
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let g2p_list = [
        G2PSchemeType::MLSMPM,
        G2PSchemeType::Lsmps4th,
        G2PSchemeType::Lsmps3rd,
        G2PSchemeType::LSMPS,
        G2PSchemeType::LsmpsLinear,
    ];

    let res_list = [50, 100, 250, 500, 1000, 2000, 4000, 8000];

    let rows = res_list.map(|res| {
        g2p_list.map(|g2p_scheme| {
            let settings = Settings {
                alpha: 0.,
                affine: true,
                space_width: 10.,
                grid_width: res,
                boundary_mirror: false,
                vx_zero: false,
                weight_type: WeightType::QuinticBSpline,
                effect_radius: 4,
                g2p_scheme: g2p_scheme,
                calc_convection_term: true,
                ..Default::default()
            };
            let mut space = new_for_taylor_green(&settings);
            space.g2p(&settings);

            let pi = std::f64::consts::PI;
            let half_domain_size = 1.;
            fn true_vel(x: f64, y: f64, u: f64, pi: f64) -> Vector2<f64> {
                Vector2::new(
                    f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
                    -f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
                )
            }

            fn true_vel_grad(x: f64, y: f64, u: f64, pi: f64) -> Matrix2<f64> {
                let k = pi / u;
                let c11 = k * f64::cos(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u);
                let c12 = -k * f64::sin(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u);
                let c21 = k * f64::sin(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u);
                let c22 = -k * f64::cos(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u);

                Matrix2::new(c11, c12, c21, c22)
            }

            fn der_2_velocity(x: f64, y: f64) -> Matrix2<f64> {
                let pi = std::f64::consts::PI;
                let u = 1.;
                Matrix2::new(
                    -pi * pi / u * f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
                    -pi * pi / u * f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
                    pi * pi / u * f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
                    pi * pi / u * f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
                )
            }

            let width = res;
            let _cell_width = settings.cell_width();
            let _grid_width = width + 1;
            let particles = space.get_particles();
            let l2_error = f64::sqrt(
                particles
                    .iter()
                    .map(|p| {
                        let x = p.x().x;
                        let y = p.x().y;
                        (p.v() - true_vel(x, y, half_domain_size, pi)).norm_squared()
                    })
                    .sum::<f64>()
                    / particles
                        .iter()
                        .map(|p| (p.x().x, p.x().y))
                        .map(|(x, y)| true_vel(x, y, half_domain_size, pi).norm_squared())
                        .sum::<f64>(),
            );

            let l1_error_grad = particles
                .iter()
                .map(|p| {
                    let x = p.x().x;
                    let y = p.x().y;
                    (p.c() - true_vel_grad(x, y, half_domain_size, pi))
                        .abs()
                        .row_sum()
                        .max()
                })
                .sum::<f64>()
                / particles
                    .iter()
                    .map(|p| (p.x().x, p.x().y))
                    .map(|(x, y)| {
                        true_vel_grad(x, y, half_domain_size, pi)
                            .abs()
                            .row_sum()
                            .max()
                    })
                    .sum::<f64>();

            let l2_error_laplacian = match g2p_scheme {
                G2PSchemeType::LSMPS | G2PSchemeType::Lsmps3rd | G2PSchemeType::Lsmps4th => {
                    f64::sqrt(
                        particles
                            .iter()
                            .map(|p| {
                                let x = p.x().x;
                                let y = p.x().y;
                                let v_lsmps = p.x_lsmps().transpose();
                                (v_lsmps
                                    .fixed_view_with_steps::<2, 2>((3, 0), (1, 0))
                                    .row_sum_tr()
                                    - der_2_velocity(x, y).column_sum())
                                .norm_squared()
                            })
                            .sum::<f64>()
                            / particles
                                .iter()
                                .map(|p| (p.x().x, p.x().y))
                                .map(|(x, y)| der_2_velocity(x, y).column_sum().norm_squared())
                                .sum::<f64>(),
                    )
                }
                _ => 0.,
            };

            (l2_error, l1_error_grad, l2_error_laplacian)
        })
    });

    {
        let file = folder.join(format!(
            "{}_vel.csv",
            current_time.format("%Y%m%d_%Hh%Mm%Ss")
        ));
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(
            &([
                vec!["res".to_string()],
                g2p_list
                    .map(|g2p_scheme| format!("{:?}", g2p_scheme))
                    .to_vec(),
            ]
            .concat()),
        )?;
        for i in 0..res_list.len() {
            writer.write_record(
                &([
                    vec![res_list[i].to_string()],
                    rows[i]
                        .map(|(l2_error, _, _)| l2_error.to_string())
                        .to_vec(),
                ]
                .concat()),
            )?;
        }
        writer.flush()?;
    }

    {
        let file = folder.join(format!(
            "{}_grad.csv",
            current_time.format("%Y%m%d_%Hh%Mm%Ss")
        ));
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(
            &([
                vec!["res".to_string()],
                g2p_list
                    .map(|g2p_scheme| format!("{:?}", g2p_scheme))
                    .to_vec(),
            ]
            .concat()),
        )?;
        for i in 0..res_list.len() {
            writer.write_record(
                &([
                    vec![res_list[i].to_string()],
                    rows[i]
                        .map(|(_, linf_error_grad, _)| linf_error_grad.to_string())
                        .to_vec(),
                ]
                .concat()),
            )?;
        }
        writer.flush()?;
    }

    {
        let file = folder.join(format!(
            "{}_laplacian.csv",
            current_time.format("%Y%m%d_%Hh%Mm%Ss")
        ));
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(
            &([
                vec!["res".to_string()],
                g2p_list
                    .map(|g2p_scheme| format!("{:?}", g2p_scheme))
                    .to_vec(),
            ]
            .concat()),
        )?;
        for i in 0..res_list.len() {
            writer.write_record(
                &([
                    vec![res_list[i].to_string()],
                    rows[i]
                        .map(|(_, _, l2_error_laplacian)| l2_error_laplacian.to_string())
                        .to_vec(),
                ]
                .concat()),
            )?;
        }
        writer.flush()?;
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
            if true {
                x += rng.gen_range(-1.0..=1.0) * p_dist * 0.3;
                y += rng.gen_range(-1.0..=1.0) * p_dist * 0.3;
            }
            let p = Particle::new(Vector2::new(x, y));
            particles.push(p);
        }
    }

    let cell_width = settings.cell_width();
    let mut grid: Vec<Node> = Vec::with_capacity((grid_width + 1) * (grid_width + 1));
    for i in 0..(grid_width + 1) * (grid_width + 1) {
        let (idx_x, idx_y) = (i % (grid_width + 1), i / (grid_width + 1));
        let (x, y) = (idx_x as f64 * cell_width, idx_y as f64 * cell_width);
        let velocity = Vector2::new(
            f64::sin(pi * (x - 5.) / half_domain_size) * f64::cos(pi * (y - 5.) / half_domain_size),
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

        grid.push(Node::new_with_vel_c((idx_x, idx_y), velocity, c));
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
