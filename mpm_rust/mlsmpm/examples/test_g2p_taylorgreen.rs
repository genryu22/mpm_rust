use mlsmpm::*;
use rand::Rng;

fn main() {
    let g2p = [
        G2PSchemeType::MLSMPM,
        G2PSchemeType::LSMPS,
        G2PSchemeType::LsmpsLinear,
    ];

    for scheme_type in g2p {
        let count = 100;
        let result = (0..count).map(|_| fun_name(P2GSchemeType::MLSMPM, scheme_type));
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
    let result = [200, 400].map(|res| {
        let settings = Settings {
            dt: 0.,
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
            p2g_scheme: p2g,
            g2p_scheme: g2p,
        };

        let mut space = new_for_taylor_green(&settings);
        space.g2p(&settings);

        let PI = std::f64::consts::PI;
        let half_domain_size = 1.;
        fn true_vel(x: f64, y: f64, U: f64, PI: f64) -> Vector2<f64> {
            Vector2::new(
                f64::sin(PI * (x - 5.) / U) * f64::cos(PI * (y - 5.) / U),
                -f64::cos(PI * (x - 5.) / U) * f64::sin(PI * (y - 5.) / U),
            )
        }

        let width = res;
        let cell_width = settings.cell_width();
        let grid_width = width + 1;
        let particles = space.get_particles();
        let l2_error = f64::sqrt(
            particles
                .iter()
                .map(|p| {
                    let x = p.x().x;
                    let y = p.x().y;
                    (p.v() - true_vel(x, y, half_domain_size, PI)).norm_squared()
                })
                .sum::<f64>()
                / particles
                    .iter()
                    .map(|p| (p.x().x, p.x().y))
                    .map(|(x, y)| true_vel(x, y, half_domain_size, PI).norm_squared())
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
            f64::sin(PI * (x - 5.) / half_domain_size) * f64::cos(PI * (y - 5.) / half_domain_size),
            -f64::cos(PI * (x - 5.) / half_domain_size)
                * f64::sin(PI * (y - 5.) / half_domain_size),
        );

        grid.push(Node::new_with_vel((idx_x, idx_y), velocity));
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
