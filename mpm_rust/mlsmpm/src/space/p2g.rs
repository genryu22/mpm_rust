use std::collections::HashMap;

use na::*;
use rayon::prelude::{IntoParallelRefMutIterator, ParallelIterator};

use crate::*;

pub fn p2g(p2g_scheme: &P2GSchemeType) -> fn(settings: &Settings, space: &mut Space) {
    match p2g_scheme {
        P2GSchemeType::MLSMPM => mlsmpm,
        P2GSchemeType::LSMPS => lsmps_2,
        P2GSchemeType::LsmpsLinear => lsmps_1,
        P2GSchemeType::Lsmps3rd => lsmps_3,
        P2GSchemeType::Lsmps4th => lsmps_4,
        P2GSchemeType::LsmpsOnlyForce => lsmps_only_force,
        P2GSchemeType::CompactLsmps => compact_2_1,
        P2GSchemeType::CompactLsmpsLinear => compact_1_1,
        P2GSchemeType::CompactOnlyVelocity => compact_only_velocity,
    }
}

fn mlsmpm(settings: &Settings, space: &mut Space) {
    for p in space.particles.iter() {
        for node in NodeMutIterator::new(
            settings,
            &mut space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
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

    for p in space.particles.iter_mut() {
        let (density, volume) = calc_density_and_volume(
            settings,
            p,
            &space.grid,
            &space.period_bounds,
            &space.period_bound_rect,
        );

        let pressure = match settings.pressure {
            Some(pressure) => pressure(p, space.steps as f64 * settings.dt),
            None => {
                let pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                    * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                if pressure < 0. {
                    0.
                } else {
                    pressure
                }
            }
        };

        p.pressure = pressure;

        let dudv = p.c;
        let strain = dudv;
        let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
        let stress = -pressure * Matrix2f::identity() + viscosity_term;
        let eq_16_term_0 = -volume
            * match settings.weight_type {
                WeightType::CubicBSpline => 3.,
                _ => 4.,
            }
            / (settings.cell_width() * settings.cell_width())
            * stress;

        for n in NodeMutIterator::new(
            settings,
            &mut space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
        ) {
            n.node.force += eq_16_term_0 * n.weight * n.dist;
        }
    }
}

mlsmpm_macro::lsmps_p2g_func!(1);
mlsmpm_macro::lsmps_p2g_func!(2);
mlsmpm_macro::lsmps_p2g_func!(3);
mlsmpm_macro::lsmps_p2g_func!(4);

fn lsmps_only_force(settings: &Settings, space: &mut Space) {
    for p in space.particles.iter() {
        for node in NodeMutIterator::new(
            settings,
            &mut space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
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

    fn poly(r: Vector2<f64>) -> Vector6<f64> {
        vector![1., r.x, r.y, r.x * r.x, r.x * r.y, r.y * r.y]
    }

    let re = settings.cell_width() * 2.;
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
        f_stress: Matrix6x3<f64>,
    }

    for p in space.particles.iter() {
        for node in NodeMutIterator::new(
            settings,
            &mut space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
        ) {
            let mass_contrib = node.weight * p.mass;
            node.node.mass += mass_contrib;
        }
    }

    let mut nodes = HashMap::new();

    for p in space.particles.iter_mut() {
        let stress = {
            let (density, volume) = calc_density_and_volume(
                settings,
                p,
                &space.grid,
                &space.period_bounds,
                &space.period_bound_rect,
            );

            let pressure = match settings.pressure {
                Some(pressure) => pressure(p, space.steps as f64 * settings.dt),
                None => {
                    let pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                        * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                    if pressure < 0. {
                        0.
                    } else {
                        pressure
                    }
                }
            };

            p.pressure = pressure;

            let dudv = p.c;
            let strain = dudv;
            let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());

            (-pressure * Matrix2f::identity() + viscosity_term)
            // volumeを含めているが、正しいか不明
        };

        {
            fn weight_function() -> fn(f64, f64) -> f64 {
                fn quadratic_b_spline(x: f64) -> f64 {
                    let x = x.abs();

                    if 0. <= x && x <= 0.5 {
                        0.75 - x * x
                    } else if 0.5 <= x && x <= 1.5 {
                        0.5 * (x - 1.5).powi(2)
                    } else {
                        0.
                    }
                }

                fn qubic_b_spline(x: f64) -> f64 {
                    let x = x.abs();

                    if 0. <= x && x <= 1. {
                        0.5 * x * x * x - x * x + 2. / 3.
                    } else if 1. <= x && x <= 2. {
                        (2. - x).powi(3) / 6.
                    } else {
                        0.
                    }
                }

                fn qubic_b_spline_2d(x: f64, y: f64) -> f64 {
                    qubic_b_spline(x) * qubic_b_spline(y)
                }

                fn quadratic_b_spline_2d(x: f64, y: f64) -> f64 {
                    quadratic_b_spline(x) * quadratic_b_spline(y)
                }

                qubic_b_spline_2d
            }

            let effect_size = 10;
            (-effect_size..=effect_size)
                .flat_map(|gx| (-effect_size..=effect_size).map(move |gy| (gx, gy)))
                .for_each(|(gx, gy)| {
                    let node_ipos =
                        p.x().map(|x| (x / settings.cell_width()).floor() as i64) + vector![gx, gy];
                    let node_pos = node_ipos.cast::<f64>() * settings.cell_width();
                    let dist = node_pos - p.x();
                    if dist.norm() > re {
                        return;
                    }
                    let node_index = if let Some(rect) = &space.period_bound_rect {
                        let x_min_index = (rect.x_min / settings.cell_width()).round() as i64;
                        let x_max_index = (rect.x_max / settings.cell_width()).round() as i64;
                        let y_min_index = (rect.y_min / settings.cell_width()).round() as i64;
                        let y_max_index = (rect.y_max / settings.cell_width()).round() as i64;

                        let origin = vector![x_min_index, y_min_index];
                        let node_ipos = node_ipos - origin;
                        let node_ipos = vector![
                            (node_ipos.x.rem_euclid(x_max_index - x_min_index) + origin.x) as U,
                            (node_ipos.y.rem_euclid(y_max_index - y_min_index) + origin.y) as U
                        ];

                        (node_ipos.x, node_ipos.y)
                    } else {
                        (node_ipos.x as U, node_ipos.y as U)
                    };

                    let params = {
                        let index = node_index;
                        if !nodes.contains_key(&index) {
                            let params = LsmpsParams {
                                m: Matrix6::<f64>::zeros(),
                                f_stress: Matrix6x3::<f64>::zeros(),
                            };
                            nodes.insert(index, params);
                        }

                        nodes.get_mut(&node_index).unwrap()
                    };
                    let r_ij = -dist / rs;
                    let poly_r_ij = poly(r_ij);
                    let weight = weight_function()(
                        dist.x / settings.cell_width(),
                        dist.y / settings.cell_width(),
                    );
                    // let weight = (1. - (dist / re).norm()).powi(2);

                    params.m += weight * poly_r_ij * poly_r_ij.transpose();
                    let stress = vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                    params.f_stress += weight * poly_r_ij.kronecker(&stress.transpose());
                });
        }
    }

    space.grid.par_iter_mut().for_each(|node| {
        if !nodes.contains_key(&node.index) {
            return;
        }
        let params = nodes.get(&node.index).unwrap();
        let m_inverse = params.m.try_inverse().unwrap();

        {
            let res = scale * m_inverse * params.f_stress;
            node.force[0] = res[(1, 0)] + res[(2, 1)];
            node.force[1] = res[(1, 1)] + res[(2, 2)];
        }
    });
}

mlsmpm_macro::compact_1_p2g_func!(1);
mlsmpm_macro::compact_1_p2g_func!(2);

fn compact_only_velocity(settings: &Settings, space: &mut Space) {
    for p in space.particles.iter() {
        for node in NodeMutIterator::new(
            settings,
            &mut space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
        ) {
            let mass_contrib = node.weight * p.mass;
            node.node.mass += mass_contrib;
        }
    }

    for p in space.particles.iter_mut() {
        let (density, volume) = calc_density_and_volume(
            settings,
            p,
            &space.grid,
            &space.period_bounds,
            &space.period_bound_rect,
        );

        let pressure = match settings.pressure {
            Some(pressure) => pressure(p, space.steps as f64 * settings.dt),
            None => {
                let pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                    * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                if pressure < 0. {
                    0.
                } else {
                    pressure
                }
            }
        };

        p.pressure = pressure;

        let dudv = p.c;
        let strain = dudv;
        let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
        let stress = -pressure * Matrix2f::identity() + viscosity_term;
        let eq_16_term_0 = -volume
            * match settings.weight_type {
                WeightType::CubicBSpline => 3.,
                _ => 4.,
            }
            / (settings.cell_width() * settings.cell_width())
            * stress;

        for n in NodeMutIterator::new(
            settings,
            &mut space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
        ) {
            n.node.force += eq_16_term_0 * n.weight * n.dist;
        }
    }

    mlsmpm_macro::lsmps_poly!(2);
    mlsmpm_macro::lsmps_params!(2);
    mlsmpm_macro::compact_lsmps_func!(2, 1);
    mlsmpm_macro::compact_lsmps_scale!(2, 1);

    let rs = settings.cell_width();
    let scale_vel = scale_2_1(rs);

    let mut nodes = HashMap::new();

    for p in space.particles.iter_mut() {
        for node in NodeIterator::new(
            settings,
            &space.grid,
            p,
            &space.period_bounds,
            &space.period_bound_rect,
        ) {
            let params = {
                let index = node.node.index;
                if !nodes.contains_key(&index) {
                    let params = LsmpsParams {
                        m: SMatrix::zeros(),
                        f_vel: SMatrix::zeros(),
                        f_stress: SMatrix::zeros(),
                        f_pressure: SMatrix::zeros(),
                    };
                    nodes.insert(index, params);
                }

                nodes.get_mut(&node.node.index).unwrap()
            };

            let r_ij = -node.dist / rs;
            let poly_r_ij = poly(r_ij);
            let weight = node.weight;

            params.m += weight * poly_r_ij * poly_r_ij.transpose();

            params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose()) * c_2_1(0, 0);
            params.f_vel += weight
                * poly_r_ij.kronecker(&p.c.column(0).transpose())
                * c_2_1(1, 0)
                * -node.dist.x;
            params.f_vel += weight
                * poly_r_ij.kronecker(&p.c.column(1).transpose())
                * c_2_1(0, 1)
                * -node.dist.y;
        }
    }

    space.grid.par_iter_mut().for_each(|node| {
        if !nodes.contains_key(&node.index) {
            return;
        }
        let params = nodes.get(&node.index).unwrap();
        if let Some(m_inverse) = (params.m + Matrix6::identity() * 0.).try_inverse() {
            let res = scale_vel * m_inverse * params.f_vel;
            node.v = res.row(0).transpose();
        }
    });
}
