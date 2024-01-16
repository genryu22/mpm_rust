/*
    粒子上にテイラーグリーン渦のt=0における速度、速度勾配の理論値を与える。
    また粘性項の理論値も求める。

    以下の２つの場合の粘性項の空間解像度に対する精度の収束性を比較する。
        1. 粒子上で計算して格子上で応力の発散を計算する
        2. 格子上で速度のラプラシアンを計算する
*/

use std::{fs, sync::Mutex};

use mlsmpm::*;
use mlsmpm_macro::*;
use nalgebra::*;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

const DYNAMIC_VISCOSITY: f64 = 1.;
const RES_LIST: [usize; 8] = [50, 100, 250, 500, 1000, 2000, 4000, 8000];

const SCHEMES: [(
    &str,
    fn(settings: &Settings) -> Vec<(f64, f64, SVector<f64, 2>)>,
); 8] = [
    ("MLSMPM", scheme_mlsmpm),
    // ("応力の発散_1_g2p(p=2)", scheme_div_force_1),
    // ("応力の発散_2_g2p(p=2)", scheme_div_force_2),
    // ("応力の発散_3_g2p(p=2)", scheme_div_force_3),
    // ("応力の発散_4_g2p(p=2)", scheme_div_force_4),
    ("応力の発散_1_g2p(p=1)", scheme_div_force_1_1st),
    ("応力の発散_2_g2p(p=1)", scheme_div_force_2_1st),
    ("応力の発散_3_g2p(p=1)", scheme_div_force_3_1st),
    ("応力の発散_4_g2p(p=1)", scheme_div_force_4_1st),
    // (
    //     "速度のラプラシアン_2_g2p(p=2)",
    //     scheme_laplacian_velocity_2_g2p_2nd,
    // ),
    // (
    //     "速度のラプラシアン_3_g2p(p=2)",
    //     scheme_laplacian_velocity_3_g2p_2nd,
    // ),
    // (
    //     "速度のラプラシアン_4_g2p(p=2)",
    //     scheme_laplacian_velocity_4_g2p_2nd,
    // ),
    ("速度のラプラシアン_2_g2p(p=1)", scheme_laplacian_velocity_2),
    ("速度のラプラシアン_3_g2p(p=1)", scheme_laplacian_velocity_3),
    ("速度のラプラシアン_4_g2p(p=1)", scheme_laplacian_velocity_4),
];

scheme_div_force!(1);
scheme_div_force!(2);
scheme_div_force!(3);
scheme_div_force!(4);
scheme_laplacian_velocity_g2p_2nd!(2);
scheme_laplacian_velocity_g2p_2nd!(3);
scheme_laplacian_velocity_g2p_2nd!(4);

scheme_div_force_1st!(1);
scheme_div_force_1st!(2);
scheme_div_force_1st!(3);
scheme_div_force_1st!(4);
scheme_laplacian_velocity!(2);
scheme_laplacian_velocity!(3);
scheme_laplacian_velocity!(4);

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let folder_name = format!("exp_taylorgreen_viscosity_simple");
    let folder = std::path::Path::new(&folder_name);
    if !folder.exists() {
        fs::create_dir(folder)?;
    }
    let current_time = chrono::Local::now();

    let result = RES_LIST.map(|res| {
        let settings = Settings {
            dynamic_viscosity: DYNAMIC_VISCOSITY,
            alpha: 0.,
            affine: true,
            space_width: 10.,
            grid_width: res,
            c: 1e1,
            boundary_mirror: false,
            vx_zero: false,
            weight_type: WeightType::QuadraticBSpline,
            effect_radius: 2,
            p2g_scheme: P2GSchemeType::MLSMPM,
            g2p_scheme: G2PSchemeType::MLSMPM,
            pressure: None,
            reset_particle_position: true,
            calc_convection_term: true,
            ..Default::default()
        };

        SCHEMES.map(|(name, scheme)| {
            let result = scheme(&settings);

            {
                let folder = folder.join("distribution");
                if !folder.exists() {
                    fs::create_dir(folder.clone()).unwrap();
                }
                let folder = folder.join(current_time.format("%Y%m%d_%Hh%Mm%Ss").to_string());
                if !folder.exists() {
                    fs::create_dir(folder.clone()).unwrap();
                }
                let folder = folder.join(name);
                if !folder.exists() {
                    fs::create_dir(folder.clone()).unwrap();
                }

                let file = folder.join(format!("{}.csv", res));
                let mut writer = csv::Writer::from_path(file).unwrap();
                writer
                    .write_record(&["x", "y", "u", "v", "u_exact", "v_exact", "l2_error_norm"])
                    .unwrap();
                for (x, y, calculated) in result.iter() {
                    let exact = viscosity_term(
                        DYNAMIC_VISCOSITY,
                        der_2_velocity(0., x.clone(), y.clone(), DYNAMIC_VISCOSITY),
                    );

                    writer
                        .write_record(&[
                            format!("{:.3}", x),
                            format!("{:.3}", y),
                            calculated.x.to_string(),
                            calculated.y.to_string(),
                            exact.x.to_string(),
                            exact.y.to_string(),
                            (calculated - exact).norm().to_string(),
                        ])
                        .unwrap();
                }
                writer.flush().unwrap();
            }

            let l2_error = f64::sqrt(
                result
                    .iter()
                    .map(|(x, y, calculated)| {
                        (calculated
                            - viscosity_term(
                                DYNAMIC_VISCOSITY,
                                der_2_velocity(0., x.clone(), y.clone(), DYNAMIC_VISCOSITY),
                            ))
                        .norm_squared()
                    })
                    .sum::<f64>()
                    / result
                        .iter()
                        .map(|(x, y, _)| {
                            viscosity_term(
                                DYNAMIC_VISCOSITY,
                                der_2_velocity(0., x.clone(), y.clone(), DYNAMIC_VISCOSITY),
                            )
                            .norm_squared()
                        })
                        .sum::<f64>(),
            );

            let linf_error = result
                .iter()
                .map(|(x, y, calculated)| {
                    (calculated
                        - viscosity_term(
                            DYNAMIC_VISCOSITY,
                            der_2_velocity(0., x.clone(), y.clone(), DYNAMIC_VISCOSITY),
                        ))
                    .norm()
                })
                .fold(0.0 / 0.0, f64::max)
                / result
                    .iter()
                    .map(|(x, y, _)| {
                        viscosity_term(
                            DYNAMIC_VISCOSITY,
                            der_2_velocity(0., x.clone(), y.clone(), DYNAMIC_VISCOSITY),
                        )
                        .norm()
                    })
                    .fold(0.0 / 0.0, f64::max);

            (l2_error, linf_error)
        })
    });

    {
        let folder = folder.join("l2");
        if !folder.exists() {
            fs::create_dir(folder.clone())?;
        }
        let file = folder.join(format!("{}.csv", current_time.format("%Y%m%d_%Hh%Mm%Ss")));
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(&([vec!["res"], SCHEMES.map(|(name, _)| name).to_vec()].concat()))?;
        for i in 0..RES_LIST.len() {
            writer.write_record(
                &([
                    vec![RES_LIST[i].to_string()],
                    result[i].map(|(l2_error, _)| l2_error.to_string()).to_vec(),
                ]
                .concat()),
            )?;
        }
        writer.flush()?;
    }

    {
        let folder = folder.join("linf");
        if !folder.exists() {
            fs::create_dir(folder.clone())?;
        }
        let file = folder.join(format!("{}.csv", current_time.format("%Y%m%d_%Hh%Mm%Ss")));
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(&([vec!["res"], SCHEMES.map(|(name, _)| name).to_vec()].concat()))?;
        for i in 0..RES_LIST.len() {
            writer.write_record(
                &([
                    vec![RES_LIST[i].to_string()],
                    result[i]
                        .map(|(_, linferror)| linferror.to_string())
                        .to_vec(),
                ]
                .concat()),
            )?;
        }
        writer.flush()?;
    }

    fn der_2_velocity(t: f64, x: f64, y: f64, nu: f64) -> Matrix2<f64> {
        let pi = std::f64::consts::PI;
        let u = 1.;
        let exp_term = f64::exp(-2. * pi * pi * t / (u * u / nu));
        Matrix2::new(
            -pi * pi / u * exp_term * f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
            -pi * pi / u * exp_term * f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
            pi * pi / u * exp_term * f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
            pi * pi / u * exp_term * f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
        )
    }

    fn viscosity_term(dynamic_viscosity: f64, der_2_v: Matrix2<f64>) -> Vector2<f64> {
        dynamic_viscosity * der_2_v.column_sum()
    }

    Ok(())
}

fn scheme_mlsmpm(settings: &Settings) -> Vec<(f64, f64, SVector<f64, 2>)> {
    let (mut particles, mut grid, periodic_boundary_rect) = new_for_taylor_green(settings);

    mlsmpm::g2p_lsmps_1st(
        &mut particles,
        &grid,
        settings,
        periodic_boundary_rect.clone(),
    );

    grid.par_iter_mut().for_each(|node| node.reset());

    let grid = grid
        .into_iter()
        .map(|node| Mutex::new(node))
        .collect::<Vec<_>>();

    let period_bounds = vec![];
    let periodic_boundary_rect = Some(periodic_boundary_rect);

    parallel!(settings, particles, |p| {
        for node in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
            let q = match settings.affine {
                true => p.c() * node.dist,
                false => SVector::zeros(),
            };
            let mass_contrib = node.weight * p.mass();
            {
                let mut node = grid[node.index].lock().unwrap();
                node.mass += mass_contrib;
                node.v += mass_contrib * (p.v() + q);
            }
        }
    });

    let force_list = {
        let mut force_list = Vec::with_capacity(grid.len());
        for _ in 0..grid.len() {
            force_list.push(Mutex::new([0.; 2]));
        }
        force_list
    };

    parallel!(settings, particles, |p| {
        let (_density, volume) =
            calc_density_and_volume(settings, p, &grid, &period_bounds, &periodic_boundary_rect);

        let dudv = p.c();
        let strain = dudv;
        let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
        let stress = viscosity_term;
        let eq_16_term_0 = -volume
            * match settings.weight_type {
                WeightType::CubicBSpline => 3.,
                _ => 4.,
            }
            / (settings.cell_width() * settings.cell_width())
            * stress;

        for n in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
            let delta = eq_16_term_0 * n.weight * n.dist;
            let mut force = force_list[n.index].lock().unwrap();
            force[0] += delta.x;
            force[1] += delta.y;
        }
    });

    let result = force_list
        .into_iter()
        .enumerate()
        .map(|(i, force)| {
            let force = force.lock().unwrap();
            Vector2::<f64>::new(force[0], force[1]) / grid[i].lock().unwrap().mass
        })
        .collect::<Vec<_>>();

    grid.iter()
        .enumerate()
        .map(|(index, _)| {
            (
                (index % (settings.grid_width + 1)) as f64 * settings.cell_width(),
                (index / (settings.grid_width + 1)) as f64 * settings.cell_width(),
                result[index],
            )
        })
        .filter(|(x, y, _)| 4. <= *x && *x < 6. && 4. <= *y && *y < 6.)
        .collect::<Vec<_>>()
}

fn new_for_taylor_green(settings: &Settings) -> (Vec<Particle>, Vec<Node>, PeriodicBoundaryRect) {
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
                x += rng.gen_range(-1.0..=1.0) * p_dist * 0.05;
                y += rng.gen_range(-1.0..=1.0) * p_dist * 0.05;
            }

            let p = Particle::new_with_mass(
                Vector2::new(x, y),
                (settings.rho_0 * (half_domain_size * 2.) * (half_domain_size * 2.))
                    / (num_x * num_x) as f64,
            );
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
        grid.push(Node::new_with_vel(
            (i % (grid_width + 1), i / (grid_width + 1)),
            velocity,
        ));
    }

    (
        particles,
        grid,
        PeriodicBoundaryRect::new(pos_x_min, pos_x_max, pos_x_min, pos_x_max),
    )
}
