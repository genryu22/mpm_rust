use crate::*;

#[derive(Debug)]
pub struct Space {
    pub(super) grid: Vec<Node>,
    pub(super) particles: Vec<Particle>,

    pub(super) slip_bounds: Vec<SlipBoundary>,
    pub(super) period_bounds: Vec<PeriodicBoundary>,
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

        Space {
            grid,
            particles,
            slip_bounds: vec![
                SlipBoundary::new(4.5, Direction::X, true),
                SlipBoundary::new(5.5, Direction::X, false),
            ],
            period_bounds: vec![PeriodicBoundary::new(
                BoundaryLine::new(4.5, true),
                BoundaryLine::new(5.5, false),
                Direction::Y,
            )],
        }
    }

    pub fn new_for_dambreak(settings: &Settings) -> Space {
        let grid_width = settings.grid_width;
        let cell_size = settings.cell_width();

        let p_dist = cell_size / 2.;

        let pos_x_min = 1.1;
        let pos_x_max = 5.1;
        let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

        let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

        for i_y in 0..num_x {
            for i_x in 0..num_x {
                let mut p = Particle::new(Vector2::new(
                    (pos_x_max - pos_x_min) * (i_x as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_x_max - pos_x_min) * (i_y as f64 + 0.5) / num_x as f64 + pos_x_min,
                ));
                p.mass = (1000. * (pos_x_max - pos_x_min) * (pos_x_max - pos_x_min))
                    / (num_x * num_x) as f64;
                particles.push(p);
            }
        }

        let mut grid: Vec<Node> = Vec::with_capacity((grid_width + 1) * (grid_width + 1));
        for _i in 0..(grid_width + 1) * (grid_width + 1) {
            grid.push(Node::new());
        }

        Space {
            grid,
            particles,
            slip_bounds: vec![
                SlipBoundary::new(1., Direction::X, true),
                SlipBoundary::new(9., Direction::X, false),
                SlipBoundary::new(1., Direction::Y, true),
                SlipBoundary::new(9., Direction::Y, false),
            ],
            period_bounds: vec![],
        }
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
                    let cell_index =
                        calc_cell_index_for_poiseuille(settings, &self.period_bounds, node_ipos);
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

            let (density, volume) = calc_density_and_volume(
                settings,
                p,
                &self.grid,
                &base_ipos,
                &weights,
                &self.period_bounds,
            );

            let mut pressure = 0.;
            if settings.c != 0. && settings.eos_power != 0. {
                pressure = settings.c / settings.eos_power
                    * ((density / 1000.).powf(settings.eos_power) - 1.);
                if pressure < 0. {
                    pressure = 0.;
                }
            }
            let pressure = pressure;

            let dudv = p.c;
            let mut strain = dudv;
            let anti_trace = strain[(1, 0)] + strain[(0, 1)];
            strain[(0, 1)] = anti_trace;
            strain[(1, 0)] = anti_trace;
            strain[(0, 0)] *= 2.;
            strain[(1, 1)] *= 2.;
            let viscosity_term = settings.dynamic_viscosity * strain;
            let stress = matrix![-pressure, 0.; 0., -pressure] + viscosity_term;
            let eq_16_term_0 =
                -volume * 4. / (settings.cell_width() * settings.cell_width()) * stress;
            for gx in 0..3 {
                for gy in 0..3 {
                    let node_ipos = base_ipos + vector![gx as U, gy as U];
                    let cell_index =
                        calc_cell_index_for_poiseuille(settings, &self.period_bounds, node_ipos);
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

            let node_pos = calc_node_pos(settings, i);
            for b in self.slip_bounds.iter() {
                let &i;

                if let Direction::X = b.direction {
                    i = &node_pos.x;
                } else {
                    i = &node_pos.y;
                }

                if b.line.calc_excess(*i) >= 0. {
                    n.v = Vector2f::zeros();
                    n.v_star = Vector2f::zeros();
                }
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
                    let cell_index =
                        calc_cell_index_for_poiseuille(settings, &self.period_bounds, node_ipos);
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

            for bound in self.period_bounds.iter() {
                let &mut i;

                if let Direction::Y = bound.direction {
                    i = &mut p.x.y;
                } else {
                    i = &mut p.x.x;
                }

                if bound.a.calc_excess(*i) > 0. {
                    *i = bound.b.plus_excess(bound.a.calc_excess(*i));
                } else if bound.b.calc_excess(*i) >= 0. {
                    *i = bound.a.plus_excess(bound.b.calc_excess(*i));
                }
            }

            for bound in self.slip_bounds.iter() {
                let &mut i;
                let &mut v_i_n; // 法線方向の速度

                if let Direction::Y = bound.direction {
                    i = &mut p.x.y;
                    v_i_n = &mut p.v.x;
                } else {
                    i = &mut p.x.x;
                    v_i_n = &mut p.v.y;
                }

                if bound.line.calc_excess(*i) >= 0. {
                    //*i = bound.line.plus_excess(bound.line.calc_excess(*i));
                    *v_i_n = 0.;
                }
            }

            p.x.x = clamp(p.x.x, 0., settings.space_width);
            p.x.y = clamp(p.x.y, 0., settings.space_width);
        }
    }
}
