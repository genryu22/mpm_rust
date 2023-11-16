use mlsmpm::*;
use rand::Rng;

fn main() {
    let p2g = [
        // P2GSchemeType::MLSMPM,
        // P2GSchemeType::LSMPS,
        // P2GSchemeType::Lsmps3rd,
        // P2GSchemeType::Lsmps4th,
        // P2GSchemeType::LsmpsLinear,
        P2GSchemeType::CompactLsmps,
        // P2GSchemeType::CompactOnlyVelocity,
        P2GSchemeType::CompactLsmpsLinear,
        P2GSchemeType::Compact_1_2,
        P2GSchemeType::Compact_2_2,
        P2GSchemeType::Compact_v_0_1,
        P2GSchemeType::Compact_v_0_2,
    ];

    for scheme_type in p2g {
        let count = 10;
        let result = (0..count).map(|_| fun_name(scheme_type, G2PSchemeType::MLSMPM));
        println!(
            "schemeType: {:?}, average = {}, (min, max) = ({}, {})",
            scheme_type,
            result.clone().sum::<f64>() / count as f64,
            result.clone().fold(0.0 / 0.0, |m, v| v.min(m)),
            result.clone().fold(0.0 / 0.0, |m, v| v.max(m))
        );
    }
}

fn fun_name(p2g: P2GSchemeType, g2p: G2PSchemeType) -> f64 {
    let result = [500, 1000].map(|res| {
        let settings = Settings {
            dt: 1e-4,
            gravity: 0.,
            dynamic_viscosity: 1e-2,
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
            g2p_scheme: g2p,
            pressure: None,
            ..Default::default()
        };

        let mut space = new_for_taylor_green(&settings);

        space.clear_grid(&settings);
        space.p2g(&settings);

        let PI = std::f64::consts::PI;
        let half_domain_size = 1.;
        fn true_vel(x: f64, y: f64, U: f64, PI: f64) -> Vector2<f64> {
            Vector2::new(
                f64::sin(PI * (x - 5.) / U) * f64::cos(PI * (y - 5.) / U),
                -f64::cos(PI * (x - 5.) / U) * f64::sin(PI * (y - 5.) / U),
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
            .filter(|(x, y, node)| {
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
                        PI,
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
                            PI,
                        )
                        .norm_squared()
                    })
                    .sum::<f64>(),
        );

        //println!("res = {}*{} l2 norm error = {}", res, res, l2_error);
        (1. / res as f64, l2_error)
    });

    //println!("{}", result[0] / result[1]);

    f64::log10(result[0].1 / result[1].1) / f64::log10(result[0].0 / result[1].0)
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
            if true {
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

            let dvxdxx = -PI * PI / (half_domain_size * half_domain_size)
                * f64::sin(PI * (x - 5.) / half_domain_size)
                * f64::cos(PI * (y - 5.) / half_domain_size);
            let dvxdxy = -PI * PI / (half_domain_size * half_domain_size)
                * f64::cos(PI * (x - 5.) / half_domain_size)
                * f64::sin(PI * (y - 5.) / half_domain_size);
            let dvxdyy = -PI * PI / (half_domain_size * half_domain_size)
                * f64::sin(PI * (x - 5.) / half_domain_size)
                * f64::cos(PI * (y - 5.) / half_domain_size);
            let dvydxx = PI * PI / (half_domain_size * half_domain_size)
                * f64::cos(PI * (x - 5.) / half_domain_size)
                * f64::sin(PI * (y - 5.) / half_domain_size);
            let dvydxy = PI * PI / (half_domain_size * half_domain_size)
                * f64::sin(PI * (x - 5.) / half_domain_size)
                * f64::cos(PI * (y - 5.) / half_domain_size);
            let dvydyy = PI * PI / (half_domain_size * half_domain_size)
                * f64::cos(PI * (x - 5.) / half_domain_size)
                * f64::sin(PI * (y - 5.) / half_domain_size);

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
