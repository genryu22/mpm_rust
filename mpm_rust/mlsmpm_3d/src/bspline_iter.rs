use nalgebra::vector;

use crate::*;

pub struct NeightborNodeMut<'a> {
    pub node: &'a mut Node,
    pub weight: f64,
    pub dist: Vector3f,
}

pub struct NodeMutIterator<'a, 'b> {
    settings: &'b Settings,
    grid: &'a mut Vec<Node>,
    base: Vector3f,
    fx: Vector3f,
    weights: [Vector3f; 3],
    gx: usize,
    gy: usize,
    gz: usize,
}

impl<'a, 'b> Iterator for NodeMutIterator<'a, 'b> {
    type Item = NeightborNodeMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.gx == 3 {
            None
        } else {
            let (weight, dist, index) = calc_weight_dist_index(
                self.settings,
                &self.base,
                &self.gx,
                &self.gy,
                &self.gz,
                &self.fx,
                &self.weights,
            );
            if let Some(index) = index {
                unsafe {
                    let node = self.grid.as_mut_ptr().add(index).as_mut().unwrap();
                    self.gy += 1;
                    if self.gy == 3 {
                        self.gy = 0;
                        self.gx += 1;
                    }
                    return Some(NeightborNodeMut { node, weight, dist });
                }
            }
            None
        }
    }
}

impl<'a, 'b, 'c, 'd> NodeMutIterator<'a, 'b> {
    pub fn new(
        settings: &'b Settings,
        grid: &'a mut Vec<Node>,
        particle: &'c Particle,
    ) -> NodeMutIterator<'a, 'b> {
        let (base, fx, weights) = calc_base_fx_weights(particle, settings);

        NodeMutIterator {
            settings,
            grid,
            base,
            fx,
            weights,
            gx: 0,
            gy: 0,
            gz: 0,
        }
    }
}
pub struct NeightborNode<'a> {
    pub node: &'a Node,
    pub weight: f64,
    pub dist: Vector3f,
}

pub struct NodeIterator<'a, 'b> {
    settings: &'b Settings,
    grid: &'a Vec<Node>,
    base: Vector3f,
    fx: Vector3f,
    weights: [Vector3f; 3],
    gx: usize,
    gy: usize,
    gz: usize,
}

impl<'a, 'b> Iterator for NodeIterator<'a, 'b> {
    type Item = NeightborNode<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.gx == 3 {
            None
        } else {
            let (weight, dist, index) = calc_weight_dist_index(
                self.settings,
                &self.base,
                &self.gx,
                &self.gy,
                &self.gz,
                &self.fx,
                &self.weights,
            );
            if let Some(index) = index {
                let node = self.grid.get(index).unwrap();
                self.gy += 1;
                if self.gy == 3 {
                    self.gy = 0;
                    self.gx += 1;
                }
                return Some(NeightborNode { node, weight, dist });
            }
            None
        }
    }
}

impl<'a, 'b, 'c, 'd> NodeIterator<'a, 'b> {
    pub fn new(
        settings: &'b Settings,
        grid: &'a Vec<Node>,
        particle: &'c Particle,
    ) -> NodeIterator<'a, 'b> {
        let (base, fx, weights) = calc_base_fx_weights(particle, settings);

        NodeIterator {
            settings,
            grid,
            base,
            fx,
            weights,
            gx: 0,
            gy: 0,
            gz: 0,
        }
    }
}

pub struct ParticleMut<'a> {
    pub particle: &'a mut Particle,
    pub density: f64,
    pub volume: f64,
}
pub struct ParticleMutIterator<'a, 'b, 'c> {
    settings: &'a Settings,
    particles: &'b mut Vec<Particle>,
    grid: &'c Vec<Node>,
    index: usize,
}

impl<'a, 'b, 'c> Iterator for ParticleMutIterator<'a, 'b, 'c> {
    type Item = ParticleMut<'b>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index == self.particles.len() {
            return None;
        }

        let particle = unsafe {
            self.particles
                .as_mut_ptr()
                .add(self.index)
                .as_mut()
                .unwrap()
        };
        self.index += 1;

        let (density, volume) = calc_density_and_volume(self.settings, particle, self.grid);

        Some(ParticleMut {
            particle,
            density,
            volume,
        })
    }
}

impl<'a, 'b, 'c> ParticleMutIterator<'a, 'b, 'c> {
    pub fn new(
        settings: &'a Settings,
        particles: &'b mut Vec<Particle>,
        grid: &'c Vec<Node>,
    ) -> ParticleMutIterator<'a, 'b, 'c> {
        Self {
            settings,
            particles,
            grid,
            index: 0,
        }
    }
}

fn calc_base_fx_weights(
    particle: &Particle,
    settings: &Settings,
) -> (Vector3f, Vector3f, [Vector3f; 3]) {
    let base = (particle.x / settings.cell_width() - vector![0.5, 0.5, 0.5]).map(|e| e.floor());
    let fx = particle.x / settings.cell_width() - base;

    let weights = {
        let w_0 = pow2(vector![1.5, 1.5, 1.5] - fx).component_mul(&Vector3f::repeat(0.5));
        let w_1 = Vector3f::repeat(0.75) - pow2(fx - vector![1., 1., 1.]);
        let w_2 = pow2(fx - vector![0.5, 0.5, 0.5]).component_mul(&Vector3f::repeat(0.5));

        [w_0, w_1, w_2]
    };

    (base, fx, weights)
}

fn calc_weight_dist_index(
    settings: &Settings,
    base: &Vector3f,
    gx: &usize,
    gy: &usize,
    gz: &usize,
    fx: &Vector3f,
    weights: &[Vector3f; 3],
) -> (f64, Vector3f, Option<usize>) {
    let node_ipos = base + vector![*gx as f64, *gy as f64, *gz as f64];
    let cell_index = calc_cell_index_for_poiseuille(settings, node_ipos);
    let weight = weights[*gx].x * weights[*gy].y;
    let dist = (vector![*gx as f64, *gy as f64, *gz as f64] - fx) * settings.cell_width();
    if cell_index <= (settings.grid_width + 1).pow(2) {
        (weight, dist, Some(cell_index))
    } else {
        (weight, dist, None)
    }
}

fn calc_cell_index_for_poiseuille(settings: &Settings, mut node_ipos: Vector3f) -> U {
    let node_ipos = node_ipos.map(|p| p.floor() as usize);

    node_ipos.x + node_ipos.y * (settings.grid_width + 1)
}

fn pow2(vec: Vector3f) -> Vector3f {
    vec.component_mul(&vec)
}

fn calc_density_and_volume(settings: &Settings, p: &Particle, grid: &Vec<Node>) -> (f64, f64) {
    let mut density = 0.;
    for n in NodeIterator::new(settings, grid, p) {
        density += n.node.mass * n.weight / (settings.cell_width() * settings.cell_width());
    }
    (density, p.mass / density)
}
