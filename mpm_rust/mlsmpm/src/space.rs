use core::num;
use std::collections::HashMap;

use na::{Matrix3, Matrix3x2, Matrix4, Matrix4x2, Matrix6, Matrix6x2, Matrix6x3, Vector3, Vector6};
use rayon::prelude::{IntoParallelRefMutIterator, ParallelIterator};

use crate::*;

mod g2p;
mod p2g;

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

    pub fn clear_grid(&mut self) {
        for n in self.grid.iter_mut() {
            n.reset();
        }
    }

    pub fn p2g(&mut self, settings: &Settings) {
        p2g::p2g(&settings.p2g_scheme)(settings, self);
    }

    pub fn update_grid(&mut self, settings: &Settings) {
        let mut vel_configs = Vec::with_capacity(self.grid.len());

        for (i, n) in self.grid.iter_mut().enumerate() {
            if n.mass <= 0. {
                continue;
            }

            match settings.p2g_scheme {
                P2GSchemeType::LSMPS
                | P2GSchemeType::LsmpsLinear
                | P2GSchemeType::Lsmps3rd
                | P2GSchemeType::Lsmps4th
                | P2GSchemeType::CompactLsmps
                | P2GSchemeType::CompactLsmpsLinear
                | P2GSchemeType::Compact_0_1
                | P2GSchemeType::Compact_0_2
                | P2GSchemeType::Compact_1_2
                | P2GSchemeType::Compact_2_2 => {
                    n.v_star = n.v
                        + settings.dt * (vector![0., settings.gravity] + n.force / settings.rho_0);
                }
                P2GSchemeType::CompactOnlyVelocity => {
                    n.v_star =
                        n.v + settings.dt * (vector![0., settings.gravity] + n.force / n.mass);
                }
                P2GSchemeType::LsmpsOnlyForce => {
                    n.v /= n.mass;
                    n.v_star = n.v
                        + settings.dt * (vector![0., settings.gravity] + n.force / settings.rho_0);
                }
                _ => {
                    n.v /= n.mass;
                    n.v_star =
                        n.v + settings.dt * (vector![0., settings.gravity] + n.force / n.mass);
                }
            };

            let node_pos = calc_node_pos(settings, i);
            for (index, b) in self.slip_bounds.iter().enumerate() {
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

                        // if index == 1 {
                        //     n.v = vector![0., -0.1];
                        //     n.v_star = vector![0., -0.1];
                        // }
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

        let current_grid = self.grid.clone();
        self.grid.par_iter_mut().for_each(|node| {
            node.c = util::calc_deriv_v(settings, node, &current_grid, &self.period_bound_rect);
        });

        if settings.reset_particle_position {
            self.particles.par_iter_mut().for_each(|p| {
                p.x = p.init_x;
            });
        }
    }

    pub fn g2p(&mut self, settings: &Settings) {
        g2p::g2p(&settings.g2p_scheme)(settings, self);

        self.particles.par_iter_mut().for_each(|p| {
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
        });
    }

    pub fn get_particle_count(&self) -> usize {
        self.particles.len()
    }

    pub fn get_nodes(&self) -> Vec<Node> {
        self.grid.clone()
    }

    pub fn get_particles(&self) -> Vec<Particle> {
        self.particles.clone()
    }

    pub fn get_steps(&self) -> usize {
        self.steps
    }
}
