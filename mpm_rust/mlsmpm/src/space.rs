use crate::*;

#[derive(Debug)]
pub struct Space {
    pub(super) grid: Vec<Node>,
    pub(super) particles: Vec<Particle>,
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
            let base_ipos = calc_base_node_ipos(settings, p.x);
            let weights = calc_weights(settings, p.x, base_ipos);
            for gx in 0..3 {
                for gy in 0..3 {
                    let node_ipos = base_ipos + vector![gx as U, gy as U];
                    let cell_index = calc_cell_index_for_poiseuille(settings, node_ipos);
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
            let base_ipos = calc_base_node_ipos(settings, p.x);
            let weights = calc_weights(settings, p.x, base_ipos);

            let (_density, volume) =
                calc_density_and_volume(settings, p, &self.grid, &base_ipos, &weights);

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
                    let cell_index = calc_cell_index_for_poiseuille(settings, node_ipos);
                    if let Some(node) = self.grid.get_mut(cell_index) {
                        let weight = weights[gx].x * weights[gy].y;
                        let node_dist = node_ipos.cast::<f64>() * settings.cell_width() - p.x;

                        node.force += eq_16_term_0 * weight * node_dist;
                    }
                }
            }
        }
    }

    pub fn update_grid(&mut self, settings: &Settings) {
        for (i, n) in self.grid.iter_mut().enumerate() {
            if n.mass <= 0. {
                continue;
            }

            n.v /= n.mass;
            n.v_star = n.v + settings.dt * (vector![0., settings.gravity] + n.force / n.mass);

            // ポアズイユ流れ boundary conditions
            let node_pos = calc_node_pos(settings, i);
            if node_pos.x <= 4.5 || node_pos.x >= 5.5 {
                n.v = Vector2f::zeros();
                n.v_star = Vector2f::zeros();
            }
        }
    }

    pub fn g2p(&mut self, settings: &Settings) {
        for p in self.particles.iter_mut() {
            let p_v_t = p.v;
            p.v = Vector2f::zeros();
            p.c = Matrix2f::zeros();

            let base_ipos = calc_base_node_ipos(settings, p.x);
            let weights = calc_weights(settings, p.x, base_ipos);
            for gx in 0..3 {
                for gy in 0..3 {
                    let node_ipos = base_ipos + vector![gx as U, gy as U];
                    let cell_index = calc_cell_index_for_poiseuille(settings, node_ipos);
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
}
