use std::vec;

use mlsmpm::*;
use nalgebra::*;
use rand::Rng;

mlsmpm_macro::lsmps_g2p_moment_func!(2);
mlsmpm_macro::lsmps_g2p_moment_func!(3);

fn main() {
    let settings = Settings {
        space_width: 10.,
        grid_width: 1000,
        weight_type: WeightType::QuarticBSpline,
        effect_radius: 5,
        ..Default::default()
    };

    {
        let moments = lsmps_2(&settings);
        let mut sum = SMatrix::zeros();
        for m in moments.iter() {
            sum = sum + m;
        }
        for row in (sum / moments.len() as f64).row_iter() {
            println!(
                "{} \\\\",
                row.iter()
                    .map(|c| format!("{:>6.3}", c))
                    .collect::<Vec<_>>()
                    .join(" & ")
            );
        }
        println!();

        let mut sum_squared = SMatrix::zeros();
        for m in moments.iter() {
            sum_squared = sum_squared + (m - sum / moments.len() as f64).map(|c| c * c);
        }
        for row in (sum_squared / moments.len() as f64).row_iter() {
            for value in row.iter() {
                print!("{:>6.3}", value);
            }
            println!();
        }
        println!();
    }

    {
        let moments = lsmps_3(&settings);
        let mut sum = SMatrix::zeros();
        for m in moments.iter() {
            sum = sum + m;
        }
        for row in (sum / moments.len() as f64).row_iter() {
            println!(
                "{} \\\\",
                row.iter()
                    .map(|c| format!("{:>6.2}", c))
                    .collect::<Vec<_>>()
                    .join(" & ")
            );
        }
        println!();

        let mut sum_squared = SMatrix::zeros();
        for m in moments.iter() {
            sum_squared = sum_squared + (m - sum / moments.len() as f64).map(|c| c * c);
        }
        for row in (sum_squared / moments.len() as f64).row_iter() {
            for value in row.iter() {
                print!("{:>6.3} ", value);
            }
            println!();
        }
        println!();
    }
}

fn new_for_taylor_green(settings: &Settings) -> (Vec<Particle>, Vec<Node>, PeriodicBoundaryRect) {
    let grid_width = settings.grid_width;

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

            let p = Particle::new(Vector2::new(x, y));
            particles.push(p);
        }
    }

    let mut grid: Vec<Node> = Vec::with_capacity((grid_width + 1) * (grid_width + 1));
    for i in 0..(grid_width + 1) * (grid_width + 1) {
        grid.push(Node::new((i % (grid_width + 1), i / (grid_width + 1))));
    }

    (
        particles,
        grid,
        PeriodicBoundaryRect::new(pos_x_min, pos_x_max, pos_x_min, pos_x_max),
    )
}
