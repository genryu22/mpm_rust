use std::vec::IntoIter;

use nalgebra::vector;

use crate::*;

pub struct NeightborNodeInfo {
    pub index: Vector3i,
    pub weight: f64,
    pub dist: Vector3f,
}

pub struct NodeIterator<'a> {
    settings: &'a Settings,
    particle_pos: Vector3f,
    indices_iter: IntoIter<(i64, i64, i64)>,
}

impl<'a> Iterator for NodeIterator<'a> {
    type Item = NeightborNodeInfo;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let next_index = self.indices_iter.next();
            if next_index.is_none() {
                return None;
            }
            let (gx, gy, gz) = next_index.unwrap();
            let target_node_index =
                self.particle_pos.map(|x| x.floor() as i64) + vector![gx, gy, gz]; //
            let node_pos = target_node_index.map(|x| x as f64) * self.settings.cell_width();
            return Some(NeightborNodeInfo {
                index: target_node_index,
                weight: calc_node_weight(&self.particle_pos, self.settings, &node_pos),
                dist: node_pos - self.particle_pos,
            });
        }
    }
}

impl<'a, 'b> NodeIterator<'a> {
    pub fn new(settings: &'a Settings, particle: &'b Particle) -> NodeIterator<'a> {
        let indices = (-3..3)
            .flat_map(|gx| {
                (-3..3).flat_map(move |gy| (-3..3).map(move |gz| (gx as i64, gy as i64, gz as i64)))
            })
            .collect::<Vec<_>>();

        NodeIterator {
            settings,
            particle_pos: particle.x,
            indices_iter: indices.into_iter(),
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
    grid: &'c Grid,
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
        grid: &'c Grid,
    ) -> ParticleMutIterator<'a, 'b, 'c> {
        Self {
            settings,
            particles,
            grid,
            index: 0,
        }
    }
}

fn quadratic_kernel(x: f64) -> f64 {
    let x = x.abs();
    if 0. <= x && x < 0.5 {
        0.75 - x * x
    } else if 0.5 <= x && 1.5 <= x {
        0.5 * (1.5 - x).powi(2)
    } else {
        0.
    }
}

fn calc_node_weight(particle_pos: &Vector3f, settings: &Settings, node_pos: &Vector3f) -> f64 {
    quadratic_kernel((particle_pos.x - node_pos.x) / settings.cell_width())
        * quadratic_kernel((particle_pos.y - node_pos.y) / settings.cell_width())
        * quadratic_kernel((particle_pos.z - node_pos.z) / settings.cell_width())
}

fn pow2(vec: Vector3f) -> Vector3f {
    vec.component_mul(&vec)
}

fn calc_density_and_volume(settings: &Settings, p: &Particle, grid: &Grid) -> (f64, f64) {
    let mut density = 0.;
    for n in NodeIterator::new(settings, p) {
        if let Some(node) = grid.get_node(n.index) {
            density += node.mass * n.weight / (settings.cell_width() * settings.cell_width());
        }
    }
    (density, p.mass / density)
}
