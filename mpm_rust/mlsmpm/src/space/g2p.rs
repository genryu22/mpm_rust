use na::*;
use rayon::prelude::{IntoParallelRefMutIterator, ParallelIterator};

use crate::*;

pub fn g2p(g2p_scheme: &G2PSchemeType) -> fn(settings: &Settings, space: &mut Space) {
    match g2p_scheme {
        G2PSchemeType::MLSMPM => mlsmpm,
        G2PSchemeType::LSMPS => lsmps_2,
        G2PSchemeType::Lsmps2ndMacro => lsmps_2,
        G2PSchemeType::LsmpsLinear => lsmps_1,
        G2PSchemeType::LsmpsLinearMacro => lsmps_1,
        G2PSchemeType::CompactLsmps => compact_lsmps,
        G2PSchemeType::Lsmps3rd => lsmps_3,
        G2PSchemeType::Lsmps4th => lsmps_4,
    }
}

fn mlsmpm(settings: &Settings, space: &mut Space) {
    parallel!(settings, space.particles, |p| {
        let p_v_t = p.v;
        p.v = Vector2f::zeros();
        p.c = Matrix2f::zeros();

        for n in NodeIndexIterator::new(settings, p, &space.period_bounds, &space.period_bound_rect)
        {
            let node = space.grid[n.index].lock().unwrap();
            p.v += (node.v_star - settings.alpha * node.v) * n.weight;

            let weighted_velocity = node.v_star * n.weight;
            p.c += weighted_velocity * n.dist.transpose();
        }

        // flip
        p.v += settings.alpha * p_v_t;

        if settings.vx_zero {
            p.v.x = 0.;
        }

        if !settings.calc_convection_term {
            p.x += p.v * settings.dt;
        }

        p.c =
            p.c * match settings.weight_type {
                WeightType::CubicBSpline => 3.,
                _ => 4.,
            } / (settings.cell_width() * settings.cell_width());
    });
}

mlsmpm_macro::lsmps_g2p_func!(1);
mlsmpm_macro::lsmps_g2p_func!(2);
mlsmpm_macro::lsmps_g2p_func!(3);
mlsmpm_macro::lsmps_g2p_func!(4);

fn compact_lsmps(settings: &Settings, space: &mut Space) {
    parallel!(settings, space.particles, |p| {
        p.v = Vector2f::zeros();
        p.c = Matrix2f::zeros();

        fn factorial(num: usize) -> f64 {
            match num {
                0 | 1 => 1.,
                _ => factorial(num - 1) * num as f64,
            }
        }

        fn c(bx: usize, by: usize, p: usize, q: usize) -> f64 {
            let b = bx + by;
            if b == 0 {
                1.
            } else if b == q {
                (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by))
                    * factorial(p)
                    / factorial(p + q)
            } else {
                (-1. as f64).powi(b as i32) * factorial(b) / (factorial(bx) * factorial(by))
                    * q as f64
                    * factorial(p + q - b)
                    / factorial(p + q)
            }
        }

        fn s(ax: usize, ay: usize, rs: f64, p: usize, q: usize) -> f64 {
            let a = ax + ay;
            let sum = [(0, 0), (1, 0), (0, 1)]
                .map(|(bx, by)| {
                    let b = bx + by;
                    if b > q || bx > ax || by > ay {
                        0.
                    } else {
                        c(bx, by, p, q) / (factorial(ax - bx) * factorial(ay - by))
                    }
                })
                .iter()
                .sum::<f64>();
            1. / (sum * rs.powi(a as i32))
        }

        fn poly(r: Vector2<f64>) -> Vector6<f64> {
            vector![1., r.x, r.y, r.x * r.x, r.x * r.y, r.y * r.y]
        }

        let _re = settings.cell_width() * 3.;
        let rs = settings.cell_width();
        let scale = Matrix6::<f64>::from_diagonal(&vector![
            s(0, 0, rs, 2, 1),
            s(1, 0, rs, 2, 1),
            s(0, 1, rs, 2, 1),
            s(2, 0, rs, 2, 1),
            s(1, 1, rs, 2, 1),
            s(0, 2, rs, 2, 1)
        ]);

        struct LsmpsParams {
            m: Matrix6<f64>,
            f_vel: Matrix6x2<f64>,
        }

        let mut params = LsmpsParams {
            m: Matrix6::<f64>::zeros(),
            f_vel: Matrix6x2::<f64>::zeros(),
        };

        for n in NodeIndexIterator::new(settings, p, &space.period_bounds, &space.period_bound_rect)
        {
            let r_ij = n.dist / rs;
            let poly_r_ij = poly(r_ij);
            let weight = n.weight;

            params.m += weight * poly_r_ij * poly_r_ij.transpose();

            {
                let node = space.grid[n.index].lock().unwrap();
                params.f_vel +=
                    weight * poly_r_ij.kronecker(&node.v_star.transpose()) * c(0, 0, 2, 1) as f64;
                params.f_vel += weight
                    * poly_r_ij.kronecker(&node.c.column(0).transpose())
                    * c(1, 0, 2, 1) as f64
                    * n.dist.x;
                params.f_vel += weight
                    * poly_r_ij.kronecker(&node.c.column(1).transpose())
                    * c(0, 1, 2, 1) as f64
                    * n.dist.y;
            }
        }

        if let Some(m_inverse) = (params.m + Matrix6::identity() * 0.).try_inverse() {
            let res = scale * m_inverse * params.f_vel;
            p.v = res.row(0).transpose();
            if settings.vx_zero {
                p.v.x = 0.;
            }
            p.x += p.v * settings.dt;
            p.c = Matrix2::new(res.m21, res.m31, res.m22, res.m32);

            p.v_lsmps = res.rows(0, res.shape().0).transpose();
        }
    });
}
