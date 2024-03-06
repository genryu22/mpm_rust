use crate::*;

pub struct NeightborNodeMut<'a> {
    pub node: &'a mut Node,
    pub weight: f64,
    pub dist: Vector2f,
}

pub struct NodeMutIterator<'a, 'b, 'c> {
    settings: &'b Settings,
    grid: &'a mut Vec<Node>,
    particle_position: Vector2f,
    period_bounds: &'c Vec<PeriodicBoundary>,
    period_bound_rect: &'c Option<PeriodicBoundaryRect>,
    radius: i32,
    gx: i32,
    gy: i32,
}

impl<'a, 'b, 'd> Iterator for NodeMutIterator<'a, 'b, 'd> {
    type Item = NeightborNodeMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.gx > self.radius {
            None
        } else {
            let (weight, dist, index) = calc_weight_dist_index(
                self.settings,
                &self.particle_position,
                &self.gx,
                &self.gy,
                self.period_bounds,
                self.period_bound_rect,
            );
            if let Some(index) = index {
                unsafe {
                    let node = self.grid.as_mut_ptr().add(index).as_mut().unwrap();
                    self.gy += 1;
                    if self.gy > self.radius {
                        self.gy = -self.radius;
                        self.gx += 1;
                    }
                    return Some(NeightborNodeMut { node, weight, dist });
                }
            }
            None
        }
    }
}

impl<'a, 'b, 'c, 'd> NodeMutIterator<'a, 'b, 'd> {
    pub fn new(
        settings: &'b Settings,
        grid: &'a mut Vec<Node>,
        particle: &'c Particle,
        period_bounds: &'d Vec<PeriodicBoundary>,
        period_bound_rect: &'d Option<PeriodicBoundaryRect>,
    ) -> NodeMutIterator<'a, 'b, 'd> {
        let radius: i32 = settings.effect_radius as i32;
        NodeMutIterator {
            settings,
            grid,
            particle_position: particle.x,
            period_bounds,
            period_bound_rect,
            radius,
            gx: -radius,
            gy: -radius,
        }
    }
}

pub struct NeightborNode<'a> {
    pub node: &'a Node,
    pub weight: f64,
    pub dist: Vector2f,
}

pub struct NodeIterator<'a, 'b, 'c> {
    settings: &'b Settings,
    grid: &'a Vec<Node>,
    particle_position: Vector2f,
    period_bounds: &'c Vec<PeriodicBoundary>,
    period_bound_rect: &'c Option<PeriodicBoundaryRect>,
    radius: i32,
    gx: i32,
    gy: i32,
}

impl<'a, 'b, 'd> Iterator for NodeIterator<'a, 'b, 'd> {
    type Item = NeightborNode<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.gx > self.radius {
            None
        } else {
            let (weight, dist, index) = calc_weight_dist_index(
                self.settings,
                &self.particle_position,
                &self.gx,
                &self.gy,
                self.period_bounds,
                self.period_bound_rect,
            );
            if let Some(index) = index {
                self.gy += 1;
                if self.gy > self.radius {
                    self.gy = -self.radius;
                    self.gx += 1;
                }
                Some(NeightborNode {
                    node: self.grid.get(index).unwrap(),
                    weight,
                    dist,
                })
            } else {
                None
            }
        }
    }
}

impl<'a, 'b, 'c, 'd> NodeIterator<'a, 'b, 'd> {
    pub fn new(
        settings: &'b Settings,
        grid: &'a Vec<Node>,
        particle: &'c Particle,
        period_bounds: &'d Vec<PeriodicBoundary>,
        period_bound_rect: &'d Option<PeriodicBoundaryRect>,
    ) -> NodeIterator<'a, 'b, 'd> {
        let radius: i32 = settings.effect_radius as i32;
        NodeIterator {
            settings,
            grid,
            particle_position: particle.x,
            period_bounds,
            period_bound_rect,
            radius,
            gx: -radius,
            gy: -radius,
        }
    }
}

pub struct NeightborNodeIndex {
    pub index: usize,
    pub weight: f64,
    pub dist: Vector2f,
}

pub struct NodeIndexIterator<'a, 'b> {
    settings: &'a Settings,
    particle_position: Vector2f,
    period_bounds: &'b Vec<PeriodicBoundary>,
    period_bound_rect: &'b Option<PeriodicBoundaryRect>,
    radius: i32,
    gx: i32,
    gy: i32,
}

impl<'a, 'b> Iterator for NodeIndexIterator<'a, 'b> {
    type Item = NeightborNodeIndex;

    fn next(&mut self) -> Option<Self::Item> {
        if self.gx > self.radius {
            None
        } else {
            let (weight, dist, index) = calc_weight_dist_index(
                self.settings,
                &self.particle_position,
                &self.gx,
                &self.gy,
                self.period_bounds,
                self.period_bound_rect,
            );
            if let Some(index) = index {
                self.gy += 1;
                if self.gy > self.radius {
                    self.gy = -self.radius;
                    self.gx += 1;
                }
                Some(NeightborNodeIndex {
                    index,
                    weight,
                    dist,
                })
            } else {
                None
            }
        }
    }
}

impl<'a, 'b, 'c> NodeIndexIterator<'a, 'b> {
    pub fn new(
        settings: &'a Settings,
        particle: &'c Particle,
        period_bounds: &'b Vec<PeriodicBoundary>,
        period_bound_rect: &'b Option<PeriodicBoundaryRect>,
    ) -> NodeIndexIterator<'a, 'b> {
        let radius: i32 = settings.effect_radius as i32;
        NodeIndexIterator {
            settings,
            particle_position: particle.x,
            period_bounds,
            period_bound_rect,
            radius,
            gx: -radius,
            gy: -radius,
        }
    }
}

pub fn calc_deriv_v(
    settings: &Settings,
    node: &Node,
    grid: &Vec<Node>,
    period_bound_rect: &Option<PeriodicBoundaryRect>,
) -> Matrix2f {
    if period_bound_rect.is_none() {
        return Matrix2f::zeros();
    }

    fn get_node(
        settings: &Settings,
        rect: &PeriodicBoundaryRect,
        node_ipos: Vector2<i64>,
    ) -> (U, U) {
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
    }

    let rect = period_bound_rect.as_ref().unwrap();

    let u = |gx, gy| {
        let index = get_node(
            settings,
            &rect,
            Vector2::new(node.index.0 as i64 + gx, node.index.1 as i64 + gy),
        );

        grid[index.0 + index.1 * (settings.grid_width + 1)].v_star
    };

    let dx = (u(1, 0) - u(-1, 0)) / (2. * settings.cell_width());
    let dy = (u(0, 1) - u(0, -1)) / (2. * settings.cell_width());

    Matrix2f::from_columns(&[dx, dy])
}

fn calc_weight_dist_index(
    settings: &Settings,
    particle_position: &Vector2f,
    gx: &i32,
    gy: &i32,
    period_bounds: &Vec<PeriodicBoundary>,
    period_bound_rect: &Option<PeriodicBoundaryRect>,
) -> (f64, Vector2f, Option<usize>) {
    let base = (particle_position / settings.cell_width()).map(|e| e.floor()); // 粒子の左下のインデックス座標
    let node_ipos = base + vector![*gx as f64, *gy as f64];
    let cell_index =
        calc_cell_index_for_poiseuille(settings, period_bounds, period_bound_rect, node_ipos);
    let dist = node_ipos * settings.cell_width() - particle_position;
    let weight = weight_function(settings)(
        settings,
        dist.x / settings.cell_width(),
        dist.y / settings.cell_width(),
    );
    if cell_index <= (settings.grid_width + 1).pow(2) {
        (weight, dist, Some(cell_index))
    } else {
        (weight, dist, None)
    }
}

fn weight_function(settings: &Settings) -> fn(settings: &Settings, f64, f64) -> f64 {
    fn quadratic_b_spline_2d(_: &Settings, x: f64, y: f64) -> f64 {
        quadratic_b_spline(x) * quadratic_b_spline(y)
    }

    fn quadratic_b_spline_2(_: &Settings, x: f64, y: f64) -> f64 {
        quadratic_b_spline_range(x, 2.) * quadratic_b_spline_range(y, 2.)
    }

    fn qubic_b_spline_2d(_: &Settings, x: f64, y: f64) -> f64 {
        qubic_b_spline(x) * qubic_b_spline(y)
    }

    fn cubic_b_spline_1_5(_: &Settings, x: f64, y: f64) -> f64 {
        cubic_b_spline_range(x, 1.5) * cubic_b_spline_range(y, 1.5)
    }

    fn quartic(_: &Settings, x: f64, y: f64) -> f64 {
        quartic_b_spline(x, 2.5) * quartic_b_spline(y, 2.5)
    }

    fn quintic(_: &Settings, x: f64, y: f64) -> f64 {
        quintic_b_spline(x, 3.) * quintic_b_spline(y, 3.)
    }

    fn hexic(_: &Settings, x: f64, y: f64) -> f64 {
        hexic_b_spline(x, 3.5) * hexic_b_spline(y, 3.5)
    }

    fn heptic(_: &Settings, x: f64, y: f64) -> f64 {
        heptic_b_spline(x, 4.) * heptic_b_spline(y, 4.)
    }

    fn linear_2d(_: &Settings, x: f64, y: f64) -> f64 {
        linear(x) * linear(y)
    }

    fn spike(settings: &Settings, x: f64, y: f64) -> f64 {
        let dist = Vector2f::new(x, y);
        let re = settings.effect_radius as f64;
        if dist.norm() > re {
            0.
        } else {
            (1. - (dist / re).norm()).powi(2)
        }
    }

    match settings.weight_type {
        WeightType::QuadraticBSpline => quadratic_b_spline_2d,
        WeightType::QuadraticBSpline2 => quadratic_b_spline_2,
        WeightType::CubicBSpline => qubic_b_spline_2d,
        WeightType::CubicBSpline1_5 => cubic_b_spline_1_5,
        WeightType::QuarticBSpline => quartic,
        WeightType::QuinticBSpline => quintic,
        WeightType::HexicBSpline => hexic,
        WeightType::HepticBSpline => heptic,
        WeightType::Linear => linear_2d,
        WeightType::Spike => spike,
    }
}

fn linear(x: f64) -> f64 {
    if -1. <= x && x <= 0. {
        1. + x
    } else if 0. <= x && x <= 1. {
        1. - x
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

fn cubic_b_spline_range(x: f64, range: f64) -> f64 {
    let knots = vec![-range, -range / 2., 0., range / 2., range];

    b_spline_basis(0, 3, &knots, x)

    // let x = x.abs() * 2. / range;

    // if 0. <= x && x <= 1. {
    //     0.5 * x * x * x - x * x + 2. / 3.
    // } else if 1. <= x && x <= 2. {
    //     (2. - x).powi(3) / 6.
    // } else {
    //     0.
    // }
}

fn heptic_b_spline(x: f64, range: f64) -> f64 {
    let knots = vec![
        -range,
        -3. * range / 4.,
        -2. * range / 4.,
        -1. * range / 4.,
        0.,
        1. * range / 4.,
        2. * range / 4.,
        3. * range / 4.,
        range,
    ];

    b_spline_basis(0, 7, &knots, x)
}

fn hexic_b_spline(x: f64, range: f64) -> f64 {
    let knots = vec![
        -range,
        -5. * range / 7.,
        -3. * range / 7.,
        -1. * range / 7.,
        1. * range / 7.,
        3. * range / 7.,
        5. * range / 7.,
        range,
    ];

    b_spline_basis(0, 6, &knots, x)
}

fn quintic_b_spline(x: f64, range: f64) -> f64 {
    let knots = vec![
        -range,
        -2. * range / 3.,
        -range / 3.,
        0.,
        range / 3.,
        2. * range / 3.,
        range,
    ];

    b_spline_basis(0, 5, &knots, x)
}

fn quartic_b_spline(x: f64, range: f64) -> f64 {
    let knots = vec![
        -range,
        -3. * range / 5.,
        -range / 5.,
        range / 5.,
        3. * range / 5.,
        range,
    ];

    b_spline_basis(0, 4, &knots, x)
}

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

fn quadratic_b_spline_range(x: f64, range: f64) -> f64 {
    let x = x.abs() * 1.5 / range;

    if 0. <= x && x <= 0.5 {
        0.75 - x * x
    } else if 0.5 <= x && x <= 1.5 {
        0.5 * (x - 1.5).powi(2)
    } else {
        0.
    }
}

fn b_spline_basis(j: usize, k: usize, knots: &Vec<f64>, t: f64) -> f64 {
    if knots.len() as i32 - k as i32 - 2 < 0 || j > knots.len() - 2 {
        panic!();
    }

    if k == 0 {
        if knots[j] <= t && t < knots[j + 1] {
            1.
        } else {
            0.
        }
    } else {
        fn w(j: usize, k: usize, knot: &Vec<f64>, t: f64) -> f64 {
            if knot[j + k] == knot[j] {
                0.
            } else {
                (t - knot[j]) / (knot[j + k] - knot[j])
            }
        }

        w(j, k, knots, t) * b_spline_basis(j, k - 1, knots, t)
            + (1. - w(j + 1, k, knots, t)) * b_spline_basis(j + 1, k - 1, knots, t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quadratic_b_spline() {
        let knots = vec![-1.5, -0.5, 0.5, 1.5];
        let quadratic = |t| b_spline_basis(0, 2, &knots, t);

        let n = 10000;
        for i in 0..n {
            let t = -2.0 + i as f64 / n as f64 * 2.0;
            // println!("f({})={}", t, quadratic(t));
            assert!(f64::abs(quadratic_b_spline(t) - quadratic(t)) < 1e-10);
        }
    }

    #[test]
    fn test_quadratic_b_spline_range() {
        let range = 3.0;
        let knots = vec![-range, -range / 3., range / 3., range];
        let quadratic = |t| b_spline_basis(0, 2, &knots, t);

        let n = 10;
        for i in 0..=n {
            let t = -range + i as f64 / n as f64 * 2. * range;
            println!("f({})={}", t, quadratic(t));
            assert!(f64::abs(quadratic_b_spline_range(t, range) - quadratic(t)) < 1e-10);
        }
    }

    #[test]
    fn test_cubic_b_spline_range() {
        let range = 3.0;
        let knots = vec![-range, -range / 2., 0., range / 2., range];
        let cubic = |t| b_spline_basis(0, 3, &knots, t);

        let n = 10;
        for i in 0..=n {
            let t = -range + i as f64 / n as f64 * 2. * range;
            println!("f({})={}", t, cubic(t));
            assert!(f64::abs(cubic_b_spline_range(t, range) - cubic(t)) < 1e-10);
        }
    }
}

pub fn calc_density_and_volume(
    settings: &Settings,
    p: &Particle,
    grid: &Vec<std::sync::Mutex<Node>>,
    period_bounds: &Vec<PeriodicBoundary>,
    period_bound_rect: &Option<PeriodicBoundaryRect>,
) -> (f64, f64) {
    let mut density = 0.;
    for n in NodeIndexIterator::new(settings, p, period_bounds, period_bound_rect) {
        density += grid[n.index].lock().unwrap().mass * n.weight
            / (settings.cell_width() * settings.cell_width());
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
    period_bound_rect: &Option<PeriodicBoundaryRect>,
    node_ipos: Vector2f,
) -> U {
    if let Some(rect) = period_bound_rect {
        let x_min_index = (rect.x_min / settings.cell_width()).round() as i64;
        let x_max_index = (rect.x_max / settings.cell_width()).round() as i64;
        let y_min_index = (rect.y_min / settings.cell_width()).round() as i64;
        let y_max_index = (rect.y_max / settings.cell_width()).round() as i64;

        let origin = vector![x_min_index, y_min_index];
        let node_ipos = node_ipos.map(|p| p.floor() as i64) - origin;
        let node_ipos = vector![
            (node_ipos.x.rem_euclid(x_max_index - x_min_index) + origin.x) as U,
            (node_ipos.y.rem_euclid(y_max_index - y_min_index) + origin.y) as U
        ];

        node_ipos.x + node_ipos.y * (settings.grid_width + 1)
    } else {
        let mut node_ipos = node_ipos.map(|p| p.floor() as usize); // 微妙？

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
}

pub fn get_opposite_node_index(
    settings: &Settings,
    index: usize,
    bound: &SlipBoundary,
) -> Option<usize> {
    let line_ipos = BoundaryLine::<i64> {
        value: (bound.line.value as f64 / settings.cell_width()).round() as i64,
        lower: bound.line.lower,
    };

    let mut node_ipos = vector![
        (index % (settings.grid_width + 1)) as i64,
        (index / (settings.grid_width + 1)) as i64
    ];

    let a;
    if let Direction::Y = bound.direction {
        a = &mut node_ipos.y;
    } else {
        a = &mut node_ipos.x;
    }

    let excess = line_ipos.calc_excess(*a);
    if excess > 0 {
        return None;
    }
    *a = line_ipos.plus_excess(excess);
    if *a < 0 || *a >= settings.grid_width as i64 + 1 {
        return None;
    }

    Some((node_ipos.x + node_ipos.y * (settings.grid_width as i64 + 1)) as usize)
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