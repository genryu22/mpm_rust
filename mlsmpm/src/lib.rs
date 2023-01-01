extern crate nalgebra as na;

use na::{clamp, vector};
pub use na::{Matrix2, Vector2};

type U = usize;

type Vector2f = Vector2<f64>;
type Vector2u = Vector2<U>;
type Matrix2f = Matrix2<f64>;

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
            self.space.clear_grid();
            self.space.distribute_mass(self.settings);
            self.space.p2g(self.settings);
            self.space.update_grid(self.settings);
            self.space.g2p(self.settings);
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
    pub grid_width: U,
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
        let cell_size = settings.cell_width();

        let p_dist = cell_size / 2.;

        let pos_x_min = 4.5;
        let pos_x_max = 5.5;
        let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

        let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

        for i_y in 0..num_x {
            for i_x in 0..num_x {
                let mut p = Particle::new(Vector2::new(
                    (pos_x_max - pos_x_min) * (i_x as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_x_max - pos_x_min) * (i_y as f64 + 0.5) / num_x as f64 + pos_x_min,
                ));
                p.mass = (1. * (pos_x_max - pos_x_min) * (pos_x_max - pos_x_min))
                    / (num_x * num_x) as f64;
                particles.push(p);
            }
        }

        let mut grid: Vec<Node> = Vec::with_capacity((grid_width + 1) * (grid_width + 1));
        for _i in 0..(grid_width + 1) * (grid_width + 1) {
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
        for p in self.particles.iter() {
            let base_ipos = Self::calc_base_node_ipos(settings, p.x);
            let weights = Self::calc_weights(settings, p.x, base_ipos);
            for gx in 0..3 {
                for gy in 0..3 {
                    let node_ipos = base_ipos + vector![gx as U, gy as U];
                    let cell_index = Self::calc_cell_index_for_poiseuille(settings, node_ipos);
                    if let Some(node) = self.grid.get_mut(cell_index) {
                        let weight = weights[gx].x * weights[gy].y;
                        let node_dist = node_ipos.cast::<f64>() * settings.cell_width() - p.x;
                        let q = match settings.affine {
                            true => p.c * node_dist,
                            false => Vector2f::zeros(),
                        };
                        let mass_contrib = weight * p.mass;
                        node.mass += mass_contrib;
                        node.v += mass_contrib * (p.v + q);
                    }
                }
            }
        }
    }

    pub fn p2g(&mut self, settings: &Settings) {
        for p in self.particles.iter() {
            let base_ipos = Self::calc_base_node_ipos(settings, p.x);
            let weights = Self::calc_weights(settings, p.x, base_ipos);

            let (_density, volume) =
                Self::calc_density_and_volume(settings, p, &self.grid, &base_ipos, &weights);

            let dudv = p.c;
            let mut strain = dudv;
            let anti_trace = strain[(1, 0)] + strain[(0, 1)];
            strain[(0, 1)] = anti_trace;
            strain[(1, 0)] = anti_trace;
            strain[(0, 0)] *= 2.;
            strain[(1, 1)] *= 2.;
            let viscosity_term = settings.dynamic_viscosity * strain;
            let stress = viscosity_term;
            let eq_16_term_0 =
                -volume * 4. / (settings.cell_width() * settings.cell_width()) * stress;
            for gx in 0..3 {
                for gy in 0..3 {
                    let node_ipos = base_ipos + vector![gx as U, gy as U];
                    let cell_index = Self::calc_cell_index_for_poiseuille(settings, node_ipos);
                    if let Some(node) = self.grid.get_mut(cell_index) {
                        let weight = weights[gx].x * weights[gy].y;
                        let node_dist = node_ipos.cast::<f64>() * settings.cell_width() - p.x;

                        node.force += eq_16_term_0 * weight * node_dist;
                    }
                }
            }
        }
    }

    fn update_grid(&mut self, settings: &Settings) {
        for (i, n) in self.grid.iter_mut().enumerate() {
            if n.mass <= 0. {
                continue;
            }

            n.v /= n.mass;
            n.v_star = n.v + settings.dt * (vector![0., settings.gravity] + n.force / n.mass);

            // ポアズイユ流れ boundary conditions
            let node_pos = Self::calc_node_pos(settings, i);
            if node_pos.x <= 4.5 || node_pos.x >= 5.5 {
                n.v = Vector2f::zeros();
                n.v_star = Vector2f::zeros();
            }
        }
    }

    fn g2p(&mut self, settings: &Settings) {
        for p in self.particles.iter_mut() {
            let p_v_t = p.v;
            p.v = Vector2f::zeros();
            p.c = Matrix2f::zeros();

            let base_ipos = Self::calc_base_node_ipos(settings, p.x);
            let weights = Self::calc_weights(settings, p.x, base_ipos);
            for gx in 0..3 {
                for gy in 0..3 {
                    let node_ipos = base_ipos + vector![gx as U, gy as U];
                    let cell_index = Self::calc_cell_index_for_poiseuille(settings, node_ipos);
                    if let Some(node) = self.grid.get_mut(cell_index) {
                        let weight = weights[gx].x * weights[gy].y;
                        let node_dist = node_ipos.cast::<f64>() * settings.cell_width() - p.x;

                        p.v += (node.v_star - settings.alpha * node.v) * weight;
                        p.x += node.v_star * weight * settings.dt;

                        let weighted_velocity = node.v_star * weight;
                        p.c += weighted_velocity * node_dist.transpose();
                    }
                }
            }

            // flip
            p.v += settings.alpha * p_v_t;

            p.c = p.c * 4. / (settings.cell_width() * settings.cell_width());

            // ポアズイユ流れ
            if p.x.y < 4.5 {
                p.x.y = 5.5 - (4.5 - p.x.y);
            }

            p.x.x = clamp(p.x.x, 0., settings.space_width);
            p.x.y = clamp(p.x.y, 0., settings.space_width);
        }
    }

    fn calc_density_and_volume(
        settings: &Settings,
        p: &Particle,
        grid: &Vec<Node>,
        base_ipos: &Vector2u,
        weights: &[Vector2f; 3],
    ) -> (f64, f64) {
        let mut density = 0.;
        for gx in 0..3 {
            for gy in 0..3 {
                let node_ipos = base_ipos + vector![gx as U, gy as U];
                let cell_index = Self::calc_cell_index_for_poiseuille(settings, node_ipos);
                if let Some(node) = grid.get(cell_index) {
                    let weight = weights[gx].x * weights[gy].y;
                    density += node.mass * weight / (settings.cell_width() * settings.cell_width());
                }
            }
        }

        (density, p.mass / density)
    }

    fn calc_base_node_ipos(settings: &Settings, x: Vector2f) -> Vector2u {
        (x / settings.cell_width()).map(|e| e.round() as U)
    }

    fn calc_weights(settings: &Settings, x: Vector2f, ipos: Vector2u) -> [Vector2f; 3] {
        let fx = x / settings.cell_width() - ipos.cast::<f64>();
        let w_0 = Self::pow2(Vector2f::repeat(1.5) - fx).component_mul(&Vector2f::repeat(0.5));
        let w_1 = Vector2f::repeat(0.75) - Self::pow2(fx - Vector2f::repeat(1.0));
        let w_2 = Self::pow2(fx - Vector2f::repeat(0.5)).component_mul(&Vector2f::repeat(0.5));
        [w_0, w_1, w_2]
    }

    fn calc_cell_index_for_poiseuille(settings: &Settings, mut node_ipos: Vector2u) -> U {
        let min_y = (4.5 / settings.cell_width()).round() as U;
        let max_y = (5.5 / settings.cell_width()).round() as U;

        if node_ipos.y < min_y {
            node_ipos.y = max_y - (min_y - node_ipos.y);
        } else if max_y <= node_ipos.y {
            node_ipos.y = min_y + (node_ipos.y - max_y);
        }

        node_ipos.x + node_ipos.y * (settings.grid_width + 1)
    }

    fn calc_node_pos(settings: &Settings, index: U) -> Vector2f {
        let x_index = index % (settings.grid_width + 1);
        let y_index = index / (settings.grid_width + 1);

        vector![
            x_index as f64 * settings.cell_width(),
            y_index as f64 * settings.cell_width()
        ]
    }

    fn pow2(vec: Vector2f) -> Vector2f {
        vec.component_mul(&vec)
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
    c: Matrix2f,
    mass: f64,
}

impl Particle {
    pub fn new(pos: Vector2f) -> Particle {
        Particle {
            x: pos,
            v: Vector2::zeros(),
            c: Matrix2f::zeros(),
            mass: 0.0,
        }
    }
}

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
