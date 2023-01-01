extern crate nalgebra as na;

use std::ops::Div;

use na::vector;
pub use na::{Matrix2, Vector2};

type Vector2f = Vector2<f64>;
type Vector2u = Vector2<u16>;

#[derive(Debug)]
pub struct Calculator<'a> {
    settings: &'a Settings,
    space: Space,
}

impl<'a> Calculator<'a> {
    pub fn new(settings: &'a Settings, space: Space) -> Calculator<'a> {
        Calculator { settings, space }
    }

    pub fn start(&mut self, count: u32) {
        for i in 0..count {
            self.space.clear_grid();
            self.space.distribute_mass(self.settings);
        }
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

impl Settings {
    pub fn cell_width(&self) -> f64 {
        self.space_width / (self.grid_width as f64)
    }
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
        let cell_size = settings.cell_width();

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

    pub fn distribute_mass(&mut self, settings: &Settings) {
        for (i, p) in self.particles.iter_mut().enumerate() {
            let base_pos = Self::calc_base_node_pos(settings, p.x);
            let weights = Self::calc_weights(settings, p.x, base_pos);
            for gx in 0..3 {
                for gy in 0..3 {
                    let weight = weights[gx][0] * weights[gy][1];
                    let node_x = (base_pos.cast::<f64>() + vector![gx as f64, gy as f64])
                        * settings.cell_width();
                    let node_dist = node_x - p.x;
                }
            }
        }
    }

    fn calc_base_node_pos(settings: &Settings, x: Vector2f) -> Vector2u {
        (x / settings.cell_width()).map(|e| e as u16)
    }

    fn calc_weights(settings: &Settings, x: Vector2f, base: Vector2u) -> [Vector2f; 3] {
        let fx = x / settings.cell_width() - base.cast::<f64>();
        let w_0 = Self::pow2(Vector2f::repeat(1.5) - fx).component_mul(&Vector2f::repeat(0.5));
        let w_1 = Vector2f::repeat(0.75) - Self::pow2(fx - Vector2f::repeat(1.0));
        let w_2 = Self::pow2(fx - Vector2f::repeat(0.5)).component_mul(&Vector2f::repeat(0.5));
        [w_0, w_1, w_2]
    }

    fn pow2(vec: Vector2f) -> Vector2f {
        vec.component_div(&vec)
    }
}

#[derive(Debug)]
pub struct Node {
    pub v: Vector2f,
    v_star: Vector2f,
    force: Vector2f,
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
    x: Vector2f,
    v: Vector2f,
    c: Matrix2<f64>,
    mass: f64,
}

impl Particle {
    pub fn new(pos: Vector2f) -> Particle {
        Particle {
            x: pos,
            v: Vector2::zeros(),
            c: Matrix2::zeros(),
            mass: 0.0,
        }
    }
}
