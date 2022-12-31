extern crate nalgebra as na;

pub use na::{Matrix2, Vector2};

#[derive(Debug)]
pub struct Calculator<'a> {
    settings: &'a Settings,
    space: Space,
}

impl<'a> Calculator<'a> {
    pub fn new(settings: &'a Settings, space: Space) -> Calculator<'a> {
        Calculator { settings, space }
    }
}

#[derive(Debug)]
pub struct Settings {
    pub dt: f64,
    pub gravity: f64,
    pub dynamic_viscosity: f64,
    pub alpha: f64,
    pub affine: bool,
    pub space_width: f64,
    pub grid_width: usize,
}

#[derive(Debug)]
pub struct Space {
    grid: Vec<Node>,
    particles: Vec<Particle>,
}

impl Space {
    pub fn new_for_poiseuille(settings: &Settings) -> Space {
        let grid_width = settings.grid_width;
        let space_width = settings.space_width;
        let cell_size = space_width / (grid_width as f64);

        let p_dist = cell_size / 2.;

        let pos_x_min = 4.5;
        let pos_x_max = 5.5;
        let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

        let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

        for i_y in 0..num_x {
            for i_x in 0..num_x {
                particles.push(Particle::new(Vector2::new(
                    (pos_x_max - pos_x_min) * (i_x as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_x_max - pos_x_min) * (i_y as f64 + 0.5) / num_x as f64 + pos_x_min,
                )));
            }
        }

        let mut grid: Vec<Node> = Vec::with_capacity((grid_width + 1) * (grid_width + 1));
        for i in 0..grid.len() {
            grid.push(Node::new());
        }

        Space { grid, particles }
    }

    pub fn clear_grid(&mut self) {
        for n in self.grid.iter_mut() {
            n.reset();
        }
    }

    pub fn distribute_mass(&mut self) {}
}

#[derive(Debug)]
pub struct Node {
    pub v: Vector2<f64>,
    v_star: Vector2<f64>,
    force: Vector2<f64>,
    mass: f64,
}

impl Node {
    pub fn new() -> Node {
        Node {
            v: Vector2::zeros(),
            v_star: Vector2::zeros(),
            force: Vector2::zeros(),
            mass: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.v = Vector2::zeros();
        self.v_star = Vector2::zeros();
        self.force = Vector2::zeros();
        self.mass = 0.0;
    }
}

#[derive(Debug)]
pub struct Particle {
    p: Vector2<f64>,
    v: Vector2<f64>,
    c: Matrix2<f64>,
    mass: f64,
}

impl Particle {
    pub fn new(pos: Vector2<f64>) -> Particle {
        Particle {
            p: pos,
            v: Vector2::zeros(),
            c: Matrix2::zeros(),
            mass: 0.0,
        }
    }
}
