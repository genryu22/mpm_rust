extern crate nalgebra as na;

use na::{clamp, vector};
pub use na::{Matrix2, Vector2};

type U = usize;

type Vector2f = Vector2<f64>;
type Vector2u = Vector2<U>;
type Matrix2f = Matrix2<f64>;

mod boundary;
pub mod file;
mod node;
mod particle;
mod settings;
mod space;
mod util;

pub use boundary::*;
pub use node::*;
pub use particle::*;
pub use settings::*;
pub use space::*;
pub use util::*;

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
            // if _i % (count / 10) == 0 {
            //     println!("{} step", _i);
            // }
        }
    }

    pub fn update(&mut self) {
        self.space.clear_grid(self.settings);
        self.space.p2g(self.settings);
        self.space.update_grid(self.settings);
        self.space.g2p(self.settings);

        self.space.steps += 1;
    }

    pub fn get_particles(&self) -> &Vec<Particle> {
        &self.space.particles
    }

    pub fn get_grid(&self) -> Vec<Node> {
        self.space.get_nodes()
    }

    pub fn get_min_velocity(&self) -> f64 {
        self.space
            .particles
            .iter()
            .map(|p| p.v.y)
            .min_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    }

    pub fn get_max_velocity(&self) -> f64 {
        self.space
            .particles
            .iter()
            .map(|p| p.v.norm())
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    }

    pub fn get_max_nodal_velocity(&self) -> f64 {
        self.space
            .get_nodes()
            .iter()
            .map(|n| n.v.norm())
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    }

    pub fn get_max_nodal_star_velocity(&self) -> f64 {
        self.space
            .get_nodes()
            .iter()
            .map(|n| n.v_star.norm())
            .max_by(|x, y| x.partial_cmp(y).unwrap())
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

#[macro_export]
macro_rules! parallel {
    ($settings:expr, $vec:expr, $func:expr) => {
        match $settings.parallel {
            true => $vec.par_iter_mut().for_each($func),
            false => $vec.iter_mut().for_each($func),
        };
    };
}

pub fn g2p_lsmps_1st(
    particles: &mut Vec<Particle>,
    grid: &Vec<Node>,
    settings: &Settings,
    periodic_boundary_rect: PeriodicBoundaryRect,
) {
    use na::*;
    use rayon::prelude::*;

    let period_bounds = vec![];
    let periodic_boundary_rect = Some(periodic_boundary_rect);
    parallel!(settings, particles, |p| {
        p.v = SVector::zeros();
        p.c = SMatrix::zeros();

        mlsmpm_macro::lsmps_poly!(1);
        mlsmpm_macro::lsmps_scale!(1);
        mlsmpm_macro::lsmps_params_g2p!(1);

        let rs = settings.cell_width();
        let scale = scale(rs);

        let mut params = LsmpsParams {
            m: SMatrix::zeros(),
            f_vel: SMatrix::zeros(),
        };

        for n in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
            let r_ij = n.dist / rs;
            let poly_r_ij = poly(r_ij);
            let weight = n.weight;

            params.m += weight * poly_r_ij * poly_r_ij.transpose();

            {
                let node = &grid[n.index];
                params.f_vel += weight * poly_r_ij.kronecker(&node.v_star.transpose());
            }
        }

        if let Some(m_inverted) = params.m.try_inverse() {
            let res = scale * m_inverted * params.f_vel;
            p.v = res.row(0).transpose();

            if settings.vx_zero {
                p.v.x = 0.;
            }
            if !settings.calc_convection_term {
                p.x += p.v * settings.dt;
            }
            p.c = res.fixed_view::<2, 2>(1, 0).transpose().into();
            p.v_lsmps = res.rows(0, res.shape().0).transpose();
        }
    });
}

pub fn g2p_lsmps_2nd(
    particles: &mut Vec<Particle>,
    grid: &Vec<Node>,
    settings: &Settings,
    periodic_boundary_rect: PeriodicBoundaryRect,
) {
    use na::*;
    use rayon::prelude::*;

    let period_bounds = vec![];
    let periodic_boundary_rect = Some(periodic_boundary_rect);
    parallel!(settings, particles, |p| {
        p.v = SVector::zeros();
        p.c = SMatrix::zeros();

        mlsmpm_macro::lsmps_poly!(2);
        mlsmpm_macro::lsmps_scale!(2);
        mlsmpm_macro::lsmps_params_g2p!(2);

        let rs = settings.cell_width();
        let scale = scale(rs);

        let mut params = LsmpsParams {
            m: SMatrix::zeros(),
            f_vel: SMatrix::zeros(),
        };

        for n in NodeIndexIterator::new(settings, p, &period_bounds, &periodic_boundary_rect) {
            let r_ij = n.dist / rs;
            let poly_r_ij = poly(r_ij);
            let weight = n.weight;

            params.m += weight * poly_r_ij * poly_r_ij.transpose();

            {
                let node = &grid[n.index];
                params.f_vel += weight * poly_r_ij.kronecker(&node.v_star.transpose());
            }
        }

        if let Some(m_inverted) = params.m.try_inverse() {
            let res = scale * m_inverted * params.f_vel;
            p.v = res.row(0).transpose();

            if settings.vx_zero {
                p.v.x = 0.;
            }
            if !settings.calc_convection_term {
                p.x += p.v * settings.dt;
            }
            p.c = res.fixed_view::<2, 2>(1, 0).transpose().into();
            p.v_lsmps = res.rows(0, res.shape().0).transpose();
        }
    });
}
