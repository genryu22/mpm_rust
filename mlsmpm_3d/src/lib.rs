pub use nalgebra::*;

type U = usize;

type Vector3f = Vector3<f64>;
type Vector3i = Vector3<i64>;
type Vector3u = Vector3<U>;
type Matrix3f = Matrix3<f64>;

mod bspline_iter;
mod grid3d;
mod node3d;
mod particle3d;
mod settings;
mod space3d;

pub use bspline_iter::*;
pub use grid3d::*;
pub use node3d::*;
pub use particle3d::*;
pub use settings::*;
pub use space3d::*;

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

        self.space.steps += 1;
    }

    pub fn get_particles(&self) -> &Vec<Particle> {
        &self.space.particles
    }

    pub fn get_grid(&self) -> Vec<Node> {
        self.space.grid.all_nodes_cloned()
    }

    pub fn get_min_velocity(&self) -> f64 {
        self.space
            .particles
            .iter()
            .map(|p| p.v.y)
            .min_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    }

    pub fn get_max_x(&self) -> f64 {
        self.space
            .particles
            .iter()
            .map(|p| p.x.x)
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    }
}
