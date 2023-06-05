use nalgebra::{clamp, vector};

use crate::*;

#[derive(Debug)]
pub struct Space {
    pub(super) grid: Vec<Node>,
    pub(super) particles: Vec<Particle>,

    pub(super) steps: usize,

    pub(super) settings: Settings,
}

impl Space {
    pub fn clear_grid(&mut self) {
        for n in self.grid.iter_mut() {
            n.reset();
        }
    }

    pub fn distribute_mass(&mut self, settings: &Settings) {
        for p in self.particles.iter() {
            for node in NodeMutIterator::new(settings, &mut self.grid, p) {
                let q = match settings.affine {
                    true => p.c * node.dist,
                    false => Vector3f::zeros(),
                };
                let mass_contrib = node.weight * p.mass;
                node.node.mass += mass_contrib;
                node.node.v += mass_contrib * (p.v + q);
            }
        }
    }

    pub fn p2g(&mut self, settings: &Settings) {
        for ParticleMut {
            particle: p,
            density,
            volume,
        } in
            ParticleMutIterator::new(settings, &mut self.particles, &self.grid).collect::<Vec<_>>()
        // イテレーション中にノードの値が変わるので先にcollectする
        {
            let mut pressure = 0.;
            if settings.c != 0. && settings.eos_power != 0. {
                pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                    * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                if pressure < 0. {
                    pressure = 0.;
                }
            }

            p.pressure = pressure;
            let pressure = pressure;

            let dudv = p.c;
            let strain = dudv;
            let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
            let stress = -pressure * Matrix3f::identity() + viscosity_term;
            let eq_16_term_0 =
                -volume * 4. / (settings.cell_width() * settings.cell_width()) * stress;

            for n in NodeMutIterator::new(settings, &mut self.grid, p) {
                n.node.force += eq_16_term_0 * n.weight * n.dist;
            }
        }
    }

    pub fn update_grid(&mut self, settings: &Settings) {
        for (i, n) in self.grid.iter_mut().enumerate() {
            if n.mass <= 0. {
                continue;
            }

            n.v /= n.mass;
            n.v_star = n.v + settings.dt * (vector![0., 0., settings.gravity] + n.force / n.mass);
        }
    }

    pub fn g2p(&mut self, settings: &Settings) {
        for p in self.particles.iter_mut() {
            let p_v_t = p.v;
            p.v = Vector3f::zeros();
            p.c = Matrix3f::zeros();

            for n in NodeIterator::new(settings, &self.grid, p) {
                p.v += (n.node.v_star - settings.alpha * n.node.v) * n.weight;
                p.x += n.node.v_star * n.weight * settings.dt;

                let weighted_velocity = n.node.v_star * n.weight;
                p.c += weighted_velocity * n.dist.transpose();
            }

            // flip
            p.v += settings.alpha * p_v_t;

            if settings.vx_zero {
                p.v.x = 0.;
            }

            p.c = p.c * 4. / (settings.cell_width() * settings.cell_width());

            let dx = settings.cell_width();
            p.x.x = clamp(p.x.x, 3. * dx, settings.space_width - 3. * dx);
            p.x.y = clamp(p.x.y, 3. * dx, settings.space_width - 3. * dx);
        }
    }
}
