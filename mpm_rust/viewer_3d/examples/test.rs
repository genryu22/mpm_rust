use mlsmpm_3d::*;
use viewer_3d::run_window_bevy;

fn main() {
    let settings = Settings {
        dt: 1e-3,
        gravity: -100.,
        dynamic_viscosity: 1e-3,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 10,
        rho_0: 1000.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: false,
        vx_zero: false,
    };

    let space = new_for_test(&settings);
    println!("{}", space.get_particle_count());
    run_window_bevy(settings.space_width, settings, space);
}

pub fn new_for_test(settings: &Settings) -> Space {
    let grid_width = settings.grid_width;
    let cell_size = settings.cell_width();

    let p_dist = cell_size / 2.;

    let pos_x_min = 0.1;
    let pos_x_max = 2.9;
    let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

    let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

    for i_z in 0..num_x {
        for i_y in 0..num_x {
            for i_x in 0..num_x {
                let mut p = Particle::new(Vector3::new(
                    (pos_x_max - pos_x_min) * (i_x as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_x_max - pos_x_min) * (i_y as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_x_max - pos_x_min) * (i_z as f64 + 0.5) / num_x as f64 + pos_x_min,
                ));
                p.mass = (settings.rho_0 * (pos_x_max - pos_x_min) * (pos_x_max - pos_x_min))
                    / (num_x * num_x) as f64;
                particles.push(p);
            }
        }
    }

    let mut nodes: Vec<Vec<Vec<Node>>> =
        vec![
            vec![vec![Node::new(NodeType::Normal); grid_width + 1]; grid_width + 1];
            grid_width + 1
        ];
    for z in 0..grid_width + 1 {
        for y in 0..grid_width + 1 {
            for x in 0..grid_width + 1 {
                if z == 0 || y == 0 || x == 0 {
                    nodes[z][y][x] = Node::new(NodeType::LeftWall);
                } else if z == grid_width || y == grid_width || x == grid_width {
                    nodes[z][y][x] = Node::new(NodeType::RightWall);
                } else if z == 1 || y == 1 || x == 1 {
                    nodes[z][y][x] = Node::new(NodeType::LeftHalf);
                } else if z == grid_width - 1 || y == grid_width - 1 || x == grid_width - 1 {
                    nodes[z][y][x] = Node::new(NodeType::RightHalf);
                }
            }
        }
    }

    Space {
        grid: Grid {
            nodes,
            slip_walls: vec![],
        },
        particles,
        steps: 0,
        settings: settings.clone(),
    }
}
