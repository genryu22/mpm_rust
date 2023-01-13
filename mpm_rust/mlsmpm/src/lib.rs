extern crate nalgebra as na;

use na::{clamp, matrix, vector};
pub use na::{Matrix2, Vector2};

type U = usize;

type Vector2f = Vector2<f64>;
type Vector2u = Vector2<U>;
type Matrix2f = Matrix2<f64>;

mod boundary;
mod node;
mod particle;
mod settings;
mod space;
mod util;

use boundary::*;
pub use node::*;
pub use particle::*;
pub use settings::*;
pub use space::*;
use util::*;

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
        for _i in 0..count {
            self.update();
            if _i % (count / 10) == 0 {
                println!("{} step", _i);
            }
        }
    }

    pub fn update(&mut self) {
        self.space.clear_grid();
        self.space.distribute_mass(self.settings);
        self.space.p2g(self.settings);
        self.space.update_grid(self.settings);
        self.space.g2p(self.settings);
    }

    pub fn get_particles(&self) -> &Vec<Particle> {
        &self.space.particles
    }

    pub fn get_grid(&self) -> &Vec<Node> {
        &self.space.grid
    }

    pub fn get_min_velocity(&self) -> f64 {
        self.space
            .particles
            .iter()
            .map(|p| p.v.y)
            .min_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    }
}
