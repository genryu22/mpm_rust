use std::collections::HashMap;

use na::*;
use rayon::prelude::*;

use crate::*;

pub fn p2g(p2g_scheme: &P2GSchemeType) -> fn(settings: &Settings, space: &mut Space) {
    match p2g_scheme {
        P2GSchemeType::MLSMPM => mlsmpm,
        P2GSchemeType::LSMPS => lsmps_2,
        P2GSchemeType::LsmpsLinear => lsmps_1,
        P2GSchemeType::Lsmps3rd => lsmps_3,
        P2GSchemeType::Lsmps4th => lsmps_4,
        P2GSchemeType::CompactLsmps => compact_2_1,
        P2GSchemeType::CompactLsmpsLinear => compact_1_1,
        P2GSchemeType::Compact_0_1 => compact_0_1,
        P2GSchemeType::Compact_3_1 => compact_3_1,
        P2GSchemeType::Compact_0_2 => compact_0_2,
        P2GSchemeType::Compact_1_2 => compact_1_2,
        P2GSchemeType::Compact_2_2 => compact_2_2,
        P2GSchemeType::Compact_3_2 => compact_3_2,
        P2GSchemeType::Compact_3_3 => compact_3_3,
        P2GSchemeType::Compact_4_3 => compact_4_3,
        P2GSchemeType::Compact_4_4 => compact_4_4,
        P2GSchemeType::Compact_v_0_1 => compact_v_0_1,
        P2GSchemeType::Compact_v_1_1 => compact_v_1_1,
        P2GSchemeType::CompactOnlyVelocity => compact_v_2_1,
        P2GSchemeType::Compact_v_3_1 => compact_3_1,
        P2GSchemeType::Compact_v_0_2 => compact_v_0_2,
        P2GSchemeType::Compact_v_1_2 => compact_v_1_2,
        P2GSchemeType::Compact_v_2_2 => compact_v_2_2,
        P2GSchemeType::Compact_v_3_2 => compact_v_3_2,
        P2GSchemeType::Compact_v_0_3 => compact_v_0_3,
        P2GSchemeType::Compact_v_1_3 => compact_v_1_3,
        P2GSchemeType::Compact_v_2_3 => compact_v_2_3,
        P2GSchemeType::Compact_v_3_3 => compact_v_3_3,
        P2GSchemeType::Compact_Laplacian_2_2 => compact_laplacian_2_2,
        P2GSchemeType::Compact_Laplacian_3_2 => compact_laplacian_3_2,
    }
}

fn mlsmpm(settings: &Settings, space: &mut Space) {
    parallel!(settings, space.particles, |p| {
        for node in
            NodeIndexIterator::new(settings, p, &space.period_bounds, &space.period_bound_rect)
        {
            let q = match settings.affine {
                true => p.c * node.dist,
                false => Vector2f::zeros(),
            };
            let mass_contrib = node.weight * p.mass;
            {
                let mut node = space.grid[node.index].lock().unwrap();
                node.mass += mass_contrib;
                node.v += mass_contrib * (p.v + q);
            }
        }
    });

    parallel!(settings, space.particles, |p| {
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
        let stress = match settings.pressure_grad {
            Some(_) => Matrix2f::zeros(),
            None => -pressure * Matrix2f::identity(),
        } + viscosity_term;
        let eq_16_term_0 = -volume
            * match settings.weight_type {
                WeightType::CubicBSpline => 3.,
                _ => 4.,
            }
            / (settings.cell_width() * settings.cell_width())
            * stress;

        for n in NodeIndexIterator::new(settings, p, &space.period_bounds, &space.period_bound_rect)
        {
            space.grid[n.index].lock().unwrap().force += eq_16_term_0 * n.weight * n.dist;
        }
    });

    if let Some(pressure_grad) = settings.pressure_grad {
        parallel!(settings, space.grid, |n| {
            let mut n = n.lock().unwrap();
            let mass = n.mass;
            let x = n.index.0 as f64 * settings.cell_width();
            let y = n.index.1 as f64 * settings.cell_width();

            n.force +=
                -pressure_grad(x, y, space.steps as f64 * settings.dt) * mass / settings.rho_0;
        });
    }
}

mlsmpm_macro::lsmps_p2g_func!(1);
mlsmpm_macro::lsmps_p2g_func!(2);
mlsmpm_macro::lsmps_p2g_func!(3);
mlsmpm_macro::lsmps_p2g_func!(4);

mlsmpm_macro::compact_p2g_func!(0, 1);
mlsmpm_macro::compact_p2g_func!(1, 1);
mlsmpm_macro::compact_p2g_func!(2, 1);
mlsmpm_macro::compact_p2g_func!(3, 1);

mlsmpm_macro::compact_p2g_func!(0, 2);
mlsmpm_macro::compact_p2g_func!(1, 2);
mlsmpm_macro::compact_p2g_func!(2, 2);
mlsmpm_macro::compact_p2g_func!(3, 2);

mlsmpm_macro::compact_p2g_func!(3, 3);
mlsmpm_macro::compact_p2g_func!(4, 3);

mlsmpm_macro::compact_p2g_func!(4, 4);

mlsmpm_macro::compact_v_p2g_func!(0, 1);
mlsmpm_macro::compact_v_p2g_func!(1, 1);
mlsmpm_macro::compact_v_p2g_func!(2, 1);
mlsmpm_macro::compact_v_p2g_func!(3, 1);

mlsmpm_macro::compact_v_p2g_func!(0, 2);
mlsmpm_macro::compact_v_p2g_func!(1, 2);
mlsmpm_macro::compact_v_p2g_func!(2, 2);
mlsmpm_macro::compact_v_p2g_func!(3, 2);

mlsmpm_macro::compact_v_p2g_func!(0, 3);
mlsmpm_macro::compact_v_p2g_func!(1, 3);
mlsmpm_macro::compact_v_p2g_func!(2, 3);
mlsmpm_macro::compact_v_p2g_func!(3, 3);

mlsmpm_macro::compact_p2g_laplacian_func!(2, 2);
mlsmpm_macro::compact_p2g_laplacian_func!(3, 2);

fn test(settings: &Settings, space: &mut Space) {
    mlsmpm_macro::lsmps_poly!(2);
    mlsmpm_macro::lsmps_scale!(2);
    mlsmpm_macro::lsmps_params!(2);

    let rs = settings.cell_width();
    let scale = scale(rs);

    let nodes = (0..space.grid.len())
        .into_par_iter()
        .map(|_| {
            std::sync::Mutex::new(LsmpsParams {
                m: SMatrix::zeros(),
                f_vel: SMatrix::zeros(),
                f_stress: SMatrix::zeros(),
                f_pressure: SVector::zeros(),
                count: 0,
            })
        })
        .collect::<Vec<_>>();

    parallel!(settings, space.particles, |p| {
        for node in
            NodeIndexIterator::new(settings, p, &space.period_bounds, &space.period_bound_rect)
        {
            let mut params = nodes[node.index].lock().unwrap();

            let r_ij = -node.dist / rs;
            let poly_r_ij = poly(r_ij);
            let weight = node.weight;

            params.m += weight * poly_r_ij * poly_r_ij.transpose();
            params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose());
            params.count += 1;
        }
    });

    parallel!(settings, space.grid, |node| {
        let mut node = node.lock().unwrap();
        let index = node.index.0 + node.index.1 * (settings.grid_width + 1);
        let params = nodes[index].lock().unwrap();

        if params.count == 0 {
            return;
        }

        if let Some(m_inverse) = (params.m).try_inverse() {
            {
                let res = scale * m_inverse * params.f_vel;
                node.v = res.row(0).transpose();

                {
                    let x = node.index.0 as f64 * settings.cell_width();
                    let y = node.index.1 as f64 * settings.cell_width();
                    node.force =
                        -settings.pressure_grad.unwrap()(x, y, space.steps as f64 * settings.dt)
                            + settings.dynamic_viscosity
                                * settings.rho_0
                                * vector![res[(3, 0)] + res[(5, 0)], res[(3, 1)] + res[(5, 1)]];
                }
            }
        } else {
            let sing_values = params.m.singular_values();
            println!("モーメント行列の逆行列が求まりませんでした。 max={} min={} cond={} (x, y) = ({} {}), steps={}",
                sing_values.max(),
                sing_values.min(),
                sing_values.max() / sing_values.min(),
                node.index.0 as f64 * settings.cell_width(),
                node.index.1 as f64 * settings.cell_width(),
                space.steps
            );
        }
    });
}
