use core::num;

use crate::*;

#[derive(Debug)]
pub struct Space {
    pub(super) grid: Vec<Node>,
    pub(super) particles: Vec<Particle>,

    pub(super) slip_bounds: Vec<SlipBoundary>,
    pub(super) period_bounds: Vec<PeriodicBoundary>,
    pub(super) period_bound_rect: Option<PeriodicBoundaryRect>,

    pub(super) steps: usize,
}

impl Space {
    pub fn new(
        grid: Vec<Node>,
        particles: Vec<Particle>,
        slip_bounds: Vec<SlipBoundary>,
        period_bounds: Vec<PeriodicBoundary>,
        period_bound_rect: Option<PeriodicBoundaryRect>,
        steps: usize,
    ) -> Self {
        Self {
            grid,
            particles,
            slip_bounds,
            period_bounds,
            period_bound_rect,
            steps,
        }
    }

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
                p.mass = (settings.rho_0 * (pos_x_max - pos_x_min) * (pos_x_max - pos_x_min))
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
                SlipBoundary::new(4.5, Direction::X, true, true, true),
                SlipBoundary::new(5.5, Direction::X, false, true, true),
            ],
            period_bounds: vec![PeriodicBoundary::new(
                BoundaryLine::new(4.5, true),
                BoundaryLine::new(5.5, false),
                Direction::Y,
            )],
            period_bound_rect: None,
            steps: 0,
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
                p.mass = (settings.rho_0 * (pos_x_max - pos_x_min) * (pos_x_max - pos_x_min))
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
                SlipBoundary::new(1., Direction::X, true, true, false),
                SlipBoundary::new(9., Direction::X, false, true, false),
                SlipBoundary::new(1., Direction::Y, true, false, false),
                SlipBoundary::new(9., Direction::Y, false, false, false),
            ],
            period_bounds: vec![],
            period_bound_rect: None,
            steps: 0,
        }
    }

    pub fn new_for_dambreak_experiment(settings: &Settings) -> Space {
        let grid_width = settings.grid_width;
        let cell_size = settings.cell_width();

        let p_dist = cell_size / 2.;

        let pos_x_min = 3.0;
        let pos_x_max = 4.0;
        let num_x = ((pos_x_max - pos_x_min) / p_dist + 0.5) as usize;

        let pos_y_min = 4.0;
        let pos_y_max = 6.0;
        let num_y = ((pos_y_max - pos_y_min) / p_dist + 0.5) as usize;

        let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

        for i_y in 0..num_y {
            for i_x in 0..num_x {
                let mut p = Particle::new(Vector2::new(
                    (pos_x_max - pos_x_min) * (i_x as f64 + 0.5) / num_x as f64 + pos_x_min,
                    (pos_y_max - pos_y_min) * (i_y as f64 + 0.5) / num_y as f64 + pos_y_min,
                ));
                p.mass = (settings.rho_0 * (pos_x_max - pos_x_min) * (pos_y_max - pos_y_min))
                    / (num_x * num_y) as f64;
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
                SlipBoundary::new(3., Direction::X, true, true, false),
                SlipBoundary::new(7., Direction::X, false, true, false),
                SlipBoundary::new(4., Direction::Y, true, false, false),
                SlipBoundary::new(6.5, Direction::Y, false, false, false),
            ],
            period_bounds: vec![],
            period_bound_rect: None,
            steps: 0,
        }
    }

    pub fn new_for_taylor_green(settings: &Settings) -> Space {
        let grid_width = settings.grid_width;

        let PI = std::f64::consts::PI;
        let half_domain_size = 1.;

        let pos_x_min = 5. - half_domain_size;
        let pos_x_max = 5. + half_domain_size;
        let num_x = (half_domain_size * 2. / (settings.cell_width() / 2.)) as usize;
        let p_dist = half_domain_size * 2. / (num_x as f64);

        let mut particles = Vec::<Particle>::with_capacity(num_x * num_x);

        for i_y in 0..num_x {
            for i_x in 0..num_x {
                let mut p = Particle::new(Vector2::new(
                    p_dist * (i_x as f64 + 0.5) as f64 + pos_x_min,
                    p_dist * (i_y as f64 + 0.5) as f64 + pos_x_min,
                ));
                p.mass = (settings.rho_0 * (half_domain_size * 2.) * (half_domain_size * 2.))
                    / (num_x * num_x) as f64;
                p.v = vector![
                    f64::sin(PI * (p.x.x - 5.) / half_domain_size)
                        * f64::cos(PI * (p.x.y - 5.) / half_domain_size),
                    -f64::cos(PI * (p.x.x - 5.) / half_domain_size)
                        * f64::sin(PI * (p.x.y - 5.) / half_domain_size)
                ];
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
            slip_bounds: vec![],
            period_bounds: vec![],
            period_bound_rect: Some(PeriodicBoundaryRect::new(
                pos_x_min, pos_x_max, pos_x_min, pos_x_max,
            )),
            steps: 0,
        }
    }

    pub fn clear_grid(&mut self) {
        for n in self.grid.iter_mut() {
            n.reset();
        }
    }

    pub fn distribute_mass(&mut self, settings: &Settings) {
        for p in self.particles.iter() {
            for node in NodeMutIterator::new(
                settings,
                &mut self.grid,
                p,
                &self.period_bounds,
                &self.period_bound_rect,
            ) {
                let q = match settings.affine {
                    true => p.c * node.dist,
                    false => Vector2f::zeros(),
                };
                let mass_contrib = node.weight * p.mass;
                node.node.mass += mass_contrib;
                node.node.v += mass_contrib * (p.v + q);
            }
        }
    }

    pub fn p2g(&mut self, settings: &Settings) {
        for p in self.particles.iter_mut() {
            let (density, volume) = calc_density_and_volume(
                settings,
                p,
                &self.grid,
                &self.period_bounds,
                &self.period_bound_rect,
            );

            let mut pressure = 0.;
            if settings.c != 0. && settings.eos_power != 0. {
                pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                    * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                if pressure < 0. {
                    pressure = 0.;
                }
            }

            // pressure = {
            //     let PI = std::f64::consts::PI;
            //     let L = 1.;
            //     let rho = 1.;
            //     let U = 1.;
            //     let nu = 1e-2;

            //     let (x, y) = (p.x.x - 5., p.x.y - 5.);

            //     rho * U * U / 4.
            //         * f64::exp(-4. * PI * PI * (self.steps as f64) * settings.dt * nu / (L * L))
            //         * (f64::cos(2. * PI * x / L) + f64::cos(2. * PI * y / L))
            // };

            p.pressure = pressure;
            let pressure = pressure;

            let dudv = p.c;
            let strain = dudv;
            let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
            let stress = -pressure * Matrix2f::identity() + viscosity_term;
            let eq_16_term_0 =
                -volume * 3. / (settings.cell_width() * settings.cell_width()) * stress;

            for n in NodeMutIterator::new(
                settings,
                &mut self.grid,
                p,
                &self.period_bounds,
                &self.period_bound_rect,
            ) {
                n.node.force += eq_16_term_0 * n.weight * n.dist;
            }
        }
    }

    pub fn update_grid(&mut self, settings: &Settings) {
        let mut vel_configs = Vec::with_capacity(self.grid.len());

        for (i, n) in self.grid.iter_mut().enumerate() {
            if n.mass <= 0. {
                continue;
            }

            n.v /= n.mass;
            n.v_star = n.v + settings.dt * (vector![0., settings.gravity] + n.force / n.mass);

            let node_pos = calc_node_pos(settings, i);
            for b in self.slip_bounds.iter() {
                if settings.boundary_mirror && b.fixed {
                    if let Some(opposite) = get_opposite_node_index(settings, i, &b) {
                        if opposite == i {
                            n.v = Vector2f::zeros();
                            n.v_star = Vector2f::zeros();
                        } else {
                            vel_configs.push((opposite, -n.v, -n.v_star, &b.direction, b.no_slip));
                        }
                    }
                } else {
                    let i;

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

        for (target, v, v_star, direction, no_slip) in vel_configs {
            if let Some(target) = self.grid.get_mut(target) {
                match direction {
                    Direction::X => {
                        if no_slip {
                            target.v = vector![0., v.y];
                            target.v_star = vector![0., v_star.y];
                        } else {
                            target.v = vector![v.x, 0.];
                            target.v_star = vector![v_star.x, 0.];
                        }
                    }
                    Direction::Y => {
                        if no_slip {
                            target.v = vector![v.x, 0.];
                            target.v_star = vector![v_star.x, 0.];
                        } else {
                            target.v = vector![0., v.y];
                            target.v_star = vector![0., v_star.y];
                        }
                    }
                }
            }
        }
    }

    pub fn g2p(&mut self, settings: &Settings) {
        for p in self.particles.iter_mut() {
            let p_v_t = p.v;
            p.v = Vector2f::zeros();
            p.c = Matrix2f::zeros();

            for n in NodeIterator::new(
                settings,
                &self.grid,
                p,
                &self.period_bounds,
                &self.period_bound_rect,
            ) {
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

            p.c = p.c * 3. / (settings.cell_width() * settings.cell_width());

            if let Some(ref rect) = self.period_bound_rect {
                if p.x.x < rect.x_min {
                    p.x.x = rect.x_max - (rect.x_min - p.x.x);
                } else if rect.x_max < p.x.x {
                    p.x.x = rect.x_min + (p.x.x - rect.x_max);
                }

                if p.x.y < rect.y_min {
                    p.x.y = rect.y_max - (rect.y_min - p.x.y);
                } else if rect.y_max < p.x.y {
                    p.x.y = rect.y_min + (p.x.y - rect.y_max);
                }
            } else {
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
                    *i = bound.line.plus_excess(bound.line.calc_excess(*i));
                    *v_i_n = 0.;
                }
            }

            let dx = settings.cell_width();
            p.x.x = clamp(p.x.x, 3. * dx, settings.space_width - 3. * dx);
            p.x.y = clamp(p.x.y, 3. * dx, settings.space_width - 3. * dx);
        }
    }

    pub fn get_particle_count(&self) -> usize {
        self.particles.len()
    }
}
