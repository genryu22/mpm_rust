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
    gx: i32,
    gy: i32,
}

impl<'a, 'b, 'd> Iterator for NodeMutIterator<'a, 'b, 'd> {
    type Item = NeightborNodeMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.gx > 3 {
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
                    if self.gy > 3 {
                        self.gy = -3;
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
        NodeMutIterator {
            settings,
            grid,
            particle_position: particle.x,
            period_bounds,
            period_bound_rect,
            gx: -3,
            gy: -3,
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
    gx: i32,
    gy: i32,
}

impl<'a, 'b, 'd> Iterator for NodeIterator<'a, 'b, 'd> {
    type Item = NeightborNode<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.gx > 3 {
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
                if self.gy > 3 {
                    self.gy = -3;
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
        NodeIterator {
            settings,
            grid,
            particle_position: particle.x,
            period_bounds,
            period_bound_rect,
            gx: -3,
            gy: -3,
        }
    }
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
        dist.x / settings.cell_width(),
        dist.y / settings.cell_width(),
    );
    if cell_index <= (settings.grid_width + 1).pow(2) {
        (weight, dist, Some(cell_index))
    } else {
        (weight, dist, None)
    }
}

fn weight_function(settings: &Settings) -> fn(f64, f64) -> f64 {
    fn quadratic_b_spline_2d(x: f64, y: f64) -> f64 {
        quadratic_b_spline(x) * quadratic_b_spline(y)
    }

    fn qubic_b_spline_2d(x: f64, y: f64) -> f64 {
        qubic_b_spline(x) * qubic_b_spline(y)
    }

    fn linear_2d(x: f64, y: f64) -> f64 {
        linear(x) * linear(y)
    }

    match settings.weight_type {
        WeightType::QuadraticBSpline => quadratic_b_spline_2d,
        WeightType::QubicBSpline => qubic_b_spline_2d,
        WeightType::Linear => linear_2d,
        _ => quadratic_b_spline_2d,
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

pub fn calc_density_and_volume(
    settings: &Settings,
    p: &Particle,
    grid: &Vec<Node>,
    period_bounds: &Vec<PeriodicBoundary>,
    period_bound_rect: &Option<PeriodicBoundaryRect>,
) -> (f64, f64) {
    let mut density = 0.;
    for n in NodeIterator::new(settings, grid, p, period_bounds, period_bound_rect) {
        density += n.node.mass * n.weight / (settings.cell_width() * settings.cell_width());
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
    mut node_ipos: Vector2f,
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
