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
