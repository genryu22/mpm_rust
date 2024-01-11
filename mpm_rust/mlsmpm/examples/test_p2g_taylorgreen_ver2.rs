use std::{error::Error, fs, path::Path};

use mlsmpm::*;
use rand::Rng;

const DYNAMIC_VISCOSITY: f64 = 1e-3;

fn main() -> Result<(), Box<dyn Error>> {
    let current_time = chrono::Local::now();
    let folder_name = "exp_p2g_taylorgreen_ver2";
    let folder = Path::new(&folder_name);
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let p2g = [
        P2GSchemeType::MLSMPM,
        P2GSchemeType::LSMPS,
        P2GSchemeType::Lsmps3rd,
        P2GSchemeType::Lsmps4th,
        P2GSchemeType::LsmpsLinear,
        P2GSchemeType::CompactLsmps,
        P2GSchemeType::CompactOnlyVelocity,
        P2GSchemeType::CompactLsmpsLinear,
        P2GSchemeType::Compact1_2,
        P2GSchemeType::Compact2_2,
        P2GSchemeType::CompactV0_1,
        P2GSchemeType::CompactV0_2,
    ];

    let res_list = [50, 100, 250, 500, 1000, 2000, 4000, 8000];

    let rows = res_list.map(|res| {
        p2g.map(|p2g| {
            let settings = Settings {
                dt: 1e-4,
                gravity: 0.,
                dynamic_viscosity: DYNAMIC_VISCOSITY,
                alpha: 0.,
                affine: true,
                space_width: 10.,
                grid_width: res,
                rho_0: 1.,
                c: 1e1,
                eos_power: 4.,
                boundary_mirror: false,
                vx_zero: false,
                weight_type: WeightType::QuadraticBSpline,
                effect_radius: 2,
                p2g_scheme: p2g,
                g2p_scheme: G2PSchemeType::MLSMPM,
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
                ..Default::default()
            };

            let mut space = new_for_taylor_green(&settings);

            space.clear_grid(&settings);
            space.p2g(&settings);

            let pi = std::f64::consts::PI;
            let half_domain_size = 1.;
            fn true_vel(x: f64, y: f64, u: f64, pi: f64) -> Vector2<f64> {
                Vector2::new(
                    f64::sin(pi * (x - 5.) / u) * f64::cos(pi * (y - 5.) / u),
                    -f64::cos(pi * (x - 5.) / u) * f64::sin(pi * (y - 5.) / u),
                )
            }

            let width = res;
            let cell_width = 10. / width as f64;
            let grid_width = width + 1;
            let nodes = space.get_nodes();
            let nodes = nodes
                .iter()
                .enumerate()
                .map(|(index, node)| (index % grid_width, index / grid_width, node))
                .filter(|(x, y, _node)| {
                    4. <= *x as f64 * cell_width
                        && *x as f64 * cell_width < 6.
                        && 4. <= *y as f64 * cell_width
                        && *y as f64 * cell_width < 6.
                })
                .collect::<Vec<_>>();
            let l2_error = f64::sqrt(
                nodes
                    .iter()
                    .map(|(x, y, node)| {
                        ((match p2g {
                            P2GSchemeType::MLSMPM => (*node).v / node.mass,
                            _ => (*node).v,
                        }) - true_vel(
                            *x as f64 * cell_width,
                            *y as f64 * cell_width,
                            half_domain_size,
                            pi,
                        ))
                        .norm_squared()
                    })
                    .sum::<f64>()
                    / nodes
                        .iter()
                        .map(|(x, y, _)| {
                            true_vel(
                                *x as f64 * cell_width,
                                *y as f64 * cell_width,
                                half_domain_size,
                                pi,
                            )
                            .norm_squared()
                        })
                        .sum::<f64>(),
            );

            l2_error
        })
    });

    {
        let file = folder.join(format!("{}.csv", current_time.format("%Y%m%d_%Hh%Mm%Ss")));
        let mut writer = csv::Writer::from_path(file)?;
        writer.write_record(
            &([
                vec!["res".to_string()],
                p2g.map(|p2g_scheme| format!("{:?}", p2g_scheme)).to_vec(),
            ]
            .concat()),
        )?;
        for i in 0..res_list.len() {
            writer.write_record(
                &([
                    vec![res_list[i].to_string()],
                    rows[i].map(|l2_error| l2_error.to_string()).to_vec(),
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

            let x_lsmps = nalgebra::Matrix2xX::from_row_slice(&[
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
