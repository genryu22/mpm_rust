use std::{error::Error, fs, path::Path};

use mlsmpm::*;
use mlsmpm_macro::lsmps_poly;
use nalgebra::Dyn;
use rand::Rng;
use rayon::prelude::*;

fn main() -> Result<(), Box<dyn Error>> {
    let time = 1e-2;

    let PI = std::f64::consts::PI;
    let half_domain_size = 1.;
    const DYNAMIC_VISCOSITY: f64 = 1e-2;
    fn true_vel(t: f64, x: f64, y: f64, U: f64, PI: f64, nu: f64) -> Vector2<f64> {
        let exp_term = f64::exp(-2. * PI * PI * t / (U * U / nu));
        Vector2::new(
            U * exp_term * f64::sin(PI * (x - 5.) / U) * f64::cos(PI * (y - 5.) / U),
            -U * exp_term * f64::cos(PI * (x - 5.) / U) * f64::sin(PI * (y - 5.) / U),
        )
    }

    const DT: f64 = 1e-4;

    const SPACE_WIDTH: f64 = 10.;
    const SMALL_WIDTH: usize = 500;
    const BIG_WIDTH: usize = 1000;

    {
        let cell_width = SPACE_WIDTH / BIG_WIDTH as f64;
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

    println!(
        "初期の粒子配置間隔を {}m から {}m に変えたときのl2エラーの収束次数",
        SPACE_WIDTH / SMALL_WIDTH as f64 / 2.,
        SPACE_WIDTH / BIG_WIDTH as f64 / 2.
    );

    [
        (P2GSchemeType::MLSMPM, G2PSchemeType::MLSMPM),
        // (P2GSchemeType::MLSMPM, G2PSchemeType::LsmpsLinear),
        // // (P2GSchemeType::LsmpsLinear, G2PSchemeType::LsmpsLinear),
        // // (
        // //     P2GSchemeType::CompactLsmpsLinear,
        // //     G2PSchemeType::LsmpsLinear,
        // // ),
        // // (P2GSchemeType::CompactLsmpsLinear, G2PSchemeType::LSMPS),
        // // (P2GSchemeType::CompactLsmpsLinear, G2PSchemeType::Lsmps3rd),
        // // // (P2GSchemeType::LSMPS, G2PSchemeType::LSMPS),
        // // // (P2GSchemeType::Lsmps3rd, G2PSchemeType::Lsmps3rd),
        // // // (P2GSchemeType::Lsmps4th, G2PSchemeType::Lsmps4th),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::LSMPS),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::CompactLsmps),
        // (P2GSchemeType::CompactLsmps, G2PSchemeType::Lsmps3rd),
        // // // (P2GSchemeType::CompactLsmps, G2PSchemeType::MLSMPM),
        // // // (P2GSchemeType::CompactLsmps, G2PSchemeType::LsmpsLinear),
        // // // (P2GSchemeType::Compact_0_1, G2PSchemeType::LsmpsLinear),
        // (P2GSchemeType::Compact_3_1, G2PSchemeType::LSMPS),
        // // // (P2GSchemeType::Compact_0_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_1_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_2_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_2_2, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_3_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_3_2, G2PSchemeType::Lsmps3rd),
        (P2GSchemeType::Compact_3_3, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::MLSMPM),
        // (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::LSMPS),
        // (P2GSchemeType::CompactOnlyVelocity, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_v_0_1, G2PSchemeType::MLSMPM),
        // (P2GSchemeType::Compact_v_0_1, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_v_3_2, G2PSchemeType::LSMPS),
        // (P2GSchemeType::Compact_v_3_2, G2PSchemeType::Lsmps3rd),
        // (P2GSchemeType::Compact_v_3_3, G2PSchemeType::Lsmps3rd),
    ]
    .par_iter()
    .map(|&(p2g_scheme, g2p_scheme)| {
        let result = [SMALL_WIDTH, BIG_WIDTH]
            .par_iter()
            .map(|&grid_width| {
                let settings = Settings {
                    dt: DT,
                    gravity: 0.,
                    dynamic_viscosity: DYNAMIC_VISCOSITY,
                    alpha: 0.,
                    affine: true,
                    space_width: SPACE_WIDTH,
                    grid_width,
                    rho_0: 1.,
                    c: 1e1,
                    eos_power: 4.,
                    boundary_mirror: false,
                    vx_zero: false,
                    weight_type: WeightType::CubicBSpline,
                    effect_radius: 2,
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
                    reset_particle_position: true,
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

                (1. / grid_width as f64, l2_error)
            })
            .collect::<Vec<_>>();

        format!(
            "{:?}_{:?}: {}",
            p2g_scheme,
            g2p_scheme,
            f64::log10(result[0].1 / result[1].1) / f64::log10(result[0].0 / result[1].0)
        )
    })
    .collect::<Vec<_>>()
    .iter()
    .for_each(|res_log| {
        println!("{}", res_log);
    });

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
