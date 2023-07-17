use std::collections::HashMap;

use na::*;
use rayon::prelude::{IntoParallelRefMutIterator, ParallelIterator};

use crate::*;

pub fn g2p(g2p_scheme: &G2PSchemeType) -> fn(settings: &Settings, space: &mut Space) {
    match g2p_scheme {
        G2PSchemeType::MLSMPM => mlsmpm,
        G2PSchemeType::LSMPS => lsmps,
        G2PSchemeType::LsmpsLinear => lsmps_linear,
    }
}

fn mlsmpm(settings: &Settings, space: &mut Space) {
    space.particles.par_iter_mut().for_each(|p| {
        let p_v_t = p.v;
        p.v = Vector2f::zeros();
        p.c = Matrix2f::zeros();

        for n in NodeIterator::new(
            settings,
            &space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
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

        p.c = p.c * 4. / (settings.cell_width() * settings.cell_width());
    });
}

fn lsmps(settings: &Settings, space: &mut Space) {
    space.particles.par_iter_mut().for_each(|p| {
        p.v = Vector2f::zeros();
        p.c = Matrix2f::zeros();

        fn poly(r: Vector2<f64>) -> Vector6<f64> {
            vector![1., r.x, r.y, r.x * r.x, r.x * r.y, r.y * r.y]
        }

        let re = settings.cell_width() * 3.;
        let rs = settings.cell_width();
        let scale = Matrix6::<f64>::from_diagonal(&vector![
            1.,
            1. / rs,
            1. / rs,
            2. / rs / rs,
            1. / rs / rs,
            2. / rs / rs
        ]);

        struct LsmpsParams {
            m: Matrix6<f64>,
            f_vel: Matrix6x2<f64>,
        }

        let mut params = LsmpsParams {
            m: Matrix6::<f64>::zeros(),
            f_vel: Matrix6x2::<f64>::zeros(),
        };

        for n in NodeIterator::new(
            settings,
            &space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
        ) {
            let r_ij = n.dist / rs;
            let poly_r_ij = poly(r_ij);
            let weight = n.weight;

            params.m += weight * poly_r_ij * poly_r_ij.transpose();
            params.f_vel += weight * poly_r_ij.kronecker(&n.node.v_star.transpose());
        }

        if let Some(m_inverted) = params.m.try_inverse() {
            let res = scale * m_inverted * params.f_vel;
            p.v = res.row(0).transpose();

            if settings.vx_zero {
                p.v.x = 0.;
            }
            p.x += res.row(0).transpose() * settings.dt;
            p.c = Matrix2::new(res.m21, res.m31, res.m22, res.m32);
            //p.c = Matrix2::new(0., 0., -0.1, 0.);
        }
    });
}

fn lsmps_linear(settings: &Settings, space: &mut Space) {
    space.particles.par_iter_mut().for_each(|p| {
        p.v = Vector2f::zeros();
        p.c = Matrix2f::zeros();

        fn poly(r: Vector2<f64>) -> Vector3<f64> {
            vector![1., r.x, r.y]
        }

        let re = settings.cell_width() * 3.;
        let rs = settings.cell_width();
        let scale = Matrix3::<f64>::from_diagonal(&vector![1., 1. / rs, 1. / rs]);

        struct LsmpsParams {
            m: Matrix3<f64>,
            f_vel: Matrix3x2<f64>,
        }

        let mut params = LsmpsParams {
            m: Matrix3::<f64>::zeros(),
            f_vel: Matrix3x2::<f64>::zeros(),
        };

        for n in NodeIterator::new(
            settings,
            &space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
        ) {
            let r_ij = n.dist / rs;
            let poly_r_ij = poly(r_ij);
            let weight = n.weight;

            params.m += weight * poly_r_ij * poly_r_ij.transpose();
            params.f_vel += weight * poly_r_ij.kronecker(&n.node.v_star.transpose());
        }

        if let Some(m_inverted) = params.m.try_inverse() {
            let res = scale * m_inverted * params.f_vel;
            //println!("{}", m_inverted);
            p.v = res.row(0).transpose();
            if settings.vx_zero {
                p.v.x = 0.;
            }
            p.x += res.row(0).transpose() * settings.dt;
            p.c = Matrix2::new(res.m21, res.m31, res.m22, res.m32);
            // println!("{}", p.c);
        }
    });
}
