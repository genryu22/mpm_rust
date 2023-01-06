use std::ops::Deref;

use crate::*;

pub fn calc_density_and_volume(
    settings: &Settings,
    p: &Particle,
    grid: &Vec<Node>,
    base_ipos: &Vector2u,
    weights: &[Vector2f; 3],
    period_bounds: &Vec<PeriodicBoundary>,
) -> (f64, f64) {
    let mut density = 0.;
    for gx in 0..3 {
        for gy in 0..3 {
            let node_ipos = base_ipos + vector![gx as U, gy as U];
            let cell_index = calc_cell_index_for_poiseuille(settings, period_bounds, node_ipos);
            if let Some(node) = grid.get(cell_index) {
                let weight = weights[gx].x * weights[gy].y;
                density += node.mass * weight / (settings.cell_width() * settings.cell_width());
            }
        }
    }

    (density, p.mass / density)
}

pub fn calc_base_node_ipos(settings: &Settings, x: Vector2f) -> Vector2u {
    (x / settings.cell_width() - vector![0.5, 0.5]).map(|e| e.floor() as U)
}

pub fn calc_weights(settings: &Settings, x: Vector2f, ipos: Vector2u) -> [Vector2f; 3] {
    let fx = x / settings.cell_width() - ipos.cast::<f64>();
    let w_0 = pow2(vector![1.5, 1.5] - fx).component_mul(&Vector2f::repeat(0.5));
    let w_1 = Vector2f::repeat(0.75) - pow2(fx - vector![1., 1.]);
    let w_2 = pow2(fx - vector![0.5, 0.5]).component_mul(&Vector2f::repeat(0.5));
    [w_0, w_1, w_2]
}

pub fn calc_cell_index_for_poiseuille(
    settings: &Settings,
    period_bounds: &Vec<PeriodicBoundary>,
    mut node_ipos: Vector2u,
) -> U {
    for boundary in period_bounds.iter() {
        let line_a_ipos = BoundaryLine::<i64> {
            value: (boundary.a.value as f64 / settings.cell_width()).round() as i64,
            lower: boundary.a.lower,
        };
        let line_b_ipos = BoundaryLine::<i64> {
            value: (boundary.b.value as f64 / settings.cell_width()).round() as i64,
            lower: boundary.b.lower,
        };

        let &mut i;
        if let Direction::Y = boundary.direction {
            i = &mut node_ipos.y;
        } else {
            i = &mut node_ipos.x;
        }

        if line_a_ipos.calc_excess(*i as i64) > 0 {
            *i = line_b_ipos.plus_excess(line_a_ipos.calc_excess(*i as i64)) as usize;
        } else if line_b_ipos.calc_excess(*i as i64) >= 0 {
            *i = line_a_ipos.plus_excess(line_b_ipos.calc_excess(*i as i64)) as usize;
        }
    }

    node_ipos.x + node_ipos.y * (settings.grid_width + 1)
}

pub fn calc_node_pos(settings: &Settings, index: U) -> Vector2f {
    let x_index = index % (settings.grid_width + 1);
    let y_index = index / (settings.grid_width + 1);

    vector![
        x_index as f64 * settings.cell_width(),
        y_index as f64 * settings.cell_width()
    ]
}

fn pow2(vec: Vector2f) -> Vector2f {
    vec.component_mul(&vec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pow2() {
        let a: Vector2f = vector![2., 5.];
        assert_eq!(pow2(a), vector![4., 25.]);

        let a: Vector2f = vector![-2., 5.];
        assert_eq!(pow2(a), vector![4., 25.]);
    }

    #[test]
    fn test_calc_node_pos() {
        let settings = Settings {
            dt: 0.05,
            gravity: 1e-2,
            dynamic_viscosity: 1e-2,
            alpha: 0.,
            affine: true,
            space_width: 10.,
            grid_width: 200,
            c: 0.,
            eos_power: 0.,
        };

        assert_eq!(calc_node_pos(&settings, 0), vector![0., 0.]);
        assert_eq!(calc_node_pos(&settings, 201 * 201 - 1), vector![10., 10.]);
        assert_eq!(calc_node_pos(&settings, 200), vector![10., 0.]);
    }

    #[test]
    fn test_calc_cell_index_for_poiseuille() {
        let settings = Settings {
            dt: 0.05,
            gravity: 1e-2,
            dynamic_viscosity: 1e-2,
            alpha: 0.,
            affine: true,
            space_width: 10.,
            grid_width: 100,
            c: 0.,
            eos_power: 0.,
        };

        let period_bounds = vec![PeriodicBoundary::new(
            BoundaryLine::new(4.5, true),
            BoundaryLine::new(5.5, false),
            Direction::Y,
        )];

        let a = vector![0, 45];
        let res_a = calc_cell_index_for_poiseuille(&settings, &period_bounds, a);
        assert_eq!(a, vector![0, 45]);
        assert_eq!(res_a, 45 * 101);

        let b = vector![0, 46];
        assert_eq!(
            calc_cell_index_for_poiseuille(&settings, &period_bounds, b),
            b.x + b.y * 101
        );

        let b = vector![0, 55];
        assert_eq!(
            calc_cell_index_for_poiseuille(&settings, &period_bounds, b),
            b.x + 45 * 101
        );

        let b = vector![0, 56];
        assert_eq!(
            calc_cell_index_for_poiseuille(&settings, &period_bounds, b),
            b.x + 46 * 101
        );
    }

    #[test]
    fn test_calc_base_node_ipos() {
        let settings = Settings {
            dt: 0.05,
            gravity: 1e-2,
            dynamic_viscosity: 1e-2,
            alpha: 0.,
            affine: true,
            space_width: 10.,
            grid_width: 100,
            c: 0.,
            eos_power: 0.,
        };

        assert_eq!(
            calc_base_node_ipos(&settings, vector![4.551, 6.]),
            vector![45, 59]
        );

        assert_eq!(
            calc_base_node_ipos(&settings, vector![4.5, 6.]),
            vector![44, 59]
        );

        assert_eq!(
            calc_base_node_ipos(&settings, vector![1.26, 3.2]),
            vector![12, 31]
        );
    }

    #[test]
    fn test_weights() {
        let settings = Settings {
            dt: 0.05,
            gravity: 1e-2,
            dynamic_viscosity: 1e-2,
            alpha: 0.,
            affine: true,
            space_width: 10.,
            grid_width: 100,
            c: 0.,
            eos_power: 0.,
        };

        let x = vector![4.6, 5.];
        let weights = calc_weights(&settings, x, calc_base_node_ipos(&settings, x));
        let mut sum = Vector2f::zeros();
        for w in weights {
            sum += w;
        }
        assert_eq!(sum, vector![1., 1.]);

        let x = vector![4.65, 5.4];
        let weights = calc_weights(&settings, x, calc_base_node_ipos(&settings, x));
        let mut sum = Vector2f::zeros();
        for w in weights {
            sum += w;
        }

        let x = vector![4.75, 5.43];
        let weights = calc_weights(&settings, x, calc_base_node_ipos(&settings, x));
        let mut sum = Vector2f::zeros();
        for w in weights {
            sum += w;
        }
        assert_eq!(sum, vector![1., 1.]);
    }
}
