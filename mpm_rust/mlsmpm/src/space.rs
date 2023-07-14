use core::num;
use std::collections::HashMap;

use na::{Matrix3, Matrix3x2, Matrix6, Matrix6x2, Matrix6x3, Vector3, Vector6};
use rayon::prelude::{IntoParallelRefMutIterator, ParallelIterator};

use crate::*;

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
        match settings.p2g_scheme {
            P2GSchemeType::MLSMPM => {
                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
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

                for p in self.particles.iter_mut() {
                    let (density, volume) = calc_density_and_volume(
                        settings,
                        p,
                        &self.grid,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    );

                    let mut pressure = 0.;
                    if settings.c != 0. && settings.eos_power != 0. {
                        pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                            * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                        if pressure < 0. {
                            pressure = 0.;
                        }
                    }

                    p.pressure = pressure;
                    let pressure = pressure;

                    let dudv = p.c;
                    let strain = dudv;
                    let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
                    let stress = -pressure * Matrix2f::identity() + viscosity_term;
                    let eq_16_term_0 =
                        -volume * 4. / (settings.cell_width() * settings.cell_width()) * stress;

                    for n in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        n.node.force += eq_16_term_0 * n.weight * n.dist;
                    }
                }
            }
            P2GSchemeType::LSMPS => {
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
                    f_vel: Matrix6x2<f64>,
                    f_stress: Matrix6x3<f64>,
                }

                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let mass_contrib = node.weight * p.mass;
                        node.node.mass += mass_contrib;
                    }
                }

                let mut nodes = HashMap::new();

                for p in self.particles.iter_mut() {
                    let stress = {
                        let (density, volume) = calc_density_and_volume(
                            settings,
                            p,
                            &self.grid,
                            &self.period_bounds,
                            &self.period_bound_rect,
                        );

                        let mut pressure = 0.;
                        if settings.c != 0. && settings.eos_power != 0. {
                            pressure = settings.rho_0 * settings.c * settings.c
                                / settings.eos_power
                                * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                            if pressure < 0. {
                                pressure = 0.;
                            }
                        }

                        p.pressure = pressure;
                        let pressure = pressure;

                        let dudv = p.c;
                        let strain = dudv;
                        let viscosity_term =
                            settings.dynamic_viscosity * (strain + strain.transpose());

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
                                    p.x().map(|x| (x / settings.cell_width()).floor() as i64)
                                        + vector![gx, gy];
                                let node_pos = node_ipos.cast::<f64>() * settings.cell_width();
                                let dist = node_pos - p.x();
                                if dist.norm() > re {
                                    return;
                                }
                                let node_index = if let Some(rect) = &self.period_bound_rect {
                                    let x_min_index =
                                        (rect.x_min / settings.cell_width()).round() as i64;
                                    let x_max_index =
                                        (rect.x_max / settings.cell_width()).round() as i64;
                                    let y_min_index =
                                        (rect.y_min / settings.cell_width()).round() as i64;
                                    let y_max_index =
                                        (rect.y_max / settings.cell_width()).round() as i64;

                                    let origin = vector![x_min_index, y_min_index];
                                    let node_ipos = node_ipos - origin;
                                    let node_ipos = vector![
                                        (node_ipos.x.rem_euclid(x_max_index - x_min_index)
                                            + origin.x)
                                            as U,
                                        (node_ipos.y.rem_euclid(y_max_index - y_min_index)
                                            + origin.y)
                                            as U
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
                                            f_vel: Matrix6x2::<f64>::zeros(),
                                            f_stress: Matrix6x3::<f64>::zeros(),
                                        };
                                        nodes.insert(index, params);
                                    }

                                    nodes.get_mut(&node_index).unwrap()
                                };
                                let r_ij = dist / rs;
                                let poly_r_ij = poly(r_ij);
                                let weight = weight_function()(
                                    dist.x / settings.cell_width(),
                                    dist.y / settings.cell_width(),
                                );
                                // let weight = (1. - (dist / re).norm()).powi(2);

                                params.m += weight * poly_r_ij * poly_r_ij.transpose();
                                params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose());
                                let stress =
                                    vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                                params.f_stress +=
                                    weight * poly_r_ij.kronecker(&stress.transpose());
                            });
                    }

                    // for node in NodeIterator::new(
                    //     settings,
                    //     &self.grid,
                    //     p,
                    //     &self.period_bounds,
                    //     &self.period_bound_rect,
                    // ) {
                    //     let params = {
                    //         let index = node.node.index;
                    //         if !nodes.contains_key(&index) {
                    //             let params = LsmpsParams {
                    //                 m: Matrix6::<f64>::zeros(),
                    //                 f_vel: Matrix6x2::<f64>::zeros(),
                    //                 f_stress: Matrix6x3::<f64>::zeros(),
                    //             };
                    //             nodes.insert(index, params);
                    //         }

                    //         nodes.get_mut(&node.node.index).unwrap()
                    //     };

                    //     let r_ij = node.dist / rs;
                    //     let poly_r_ij = poly(r_ij);
                    //     let weight = node.weight;

                    //     params.m += weight * poly_r_ij * poly_r_ij.transpose();
                    //     params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose());
                    //     let stress = vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                    //     params.f_stress += weight * poly_r_ij.kronecker(&stress.transpose());
                    // }
                }

                self.grid.par_iter_mut().for_each(|node| {
                    if !nodes.contains_key(&node.index) {
                        return;
                    }
                    let params = nodes.get(&node.index).unwrap();
                    let m_inverse = params.m.try_inverse().unwrap();

                    {
                        let res = scale * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                    }

                    {
                        let res = scale * m_inverse * params.f_stress;
                        node.force[0] = res[(1, 0)] + res[(2, 1)];
                        node.force[1] = res[(1, 1)] + res[(2, 2)];
                    }
                });
            }
            P2GSchemeType::LsmpsLinear => {
                fn poly(r: Vector2<f64>) -> Vector3<f64> {
                    vector![1., r.x, r.y]
                }

                let re = settings.cell_width() * 3.;
                let rs = settings.cell_width();
                let scale = Matrix3::<f64>::from_diagonal(&vector![1., 1. / rs, 1. / rs]);

                struct LsmpsParams {
                    m: Matrix3<f64>,
                    f_vel: Matrix3x2<f64>,
                    f_stress: Matrix3<f64>,
                }

                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let mass_contrib = node.weight * p.mass;
                        node.node.mass += mass_contrib;
                    }
                }

                let mut nodes = HashMap::new();

                for p in self.particles.iter_mut() {
                    let stress = {
                        let (density, volume) = calc_density_and_volume(
                            settings,
                            p,
                            &self.grid,
                            &self.period_bounds,
                            &self.period_bound_rect,
                        );

                        let mut pressure = 0.;
                        if settings.c != 0. && settings.eos_power != 0. {
                            pressure = settings.rho_0 * settings.c * settings.c
                                / settings.eos_power
                                * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                            if pressure < 0. {
                                pressure = 0.;
                            }
                        }

                        p.pressure = pressure;
                        let pressure = pressure;

                        let dudv = p.c;
                        let strain = dudv;
                        let viscosity_term =
                            settings.dynamic_viscosity * (strain + strain.transpose());

                        (-pressure * Matrix2f::identity() + viscosity_term)
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
                            fn quadratic_b_spline_2d(x: f64, y: f64) -> f64 {
                                quadratic_b_spline(x) * quadratic_b_spline(y)
                            }

                            quadratic_b_spline_2d
                        }

                        let effect_size = 3;
                        (-effect_size..=effect_size)
                            .flat_map(|gx| (-effect_size..=effect_size).map(move |gy| (gx, gy)))
                            .for_each(|(gx, gy)| {
                                let node_ipos =
                                    p.x().map(|x| (x / settings.cell_width()).floor() as i64)
                                        + vector![gx, gy];
                                let node_pos = node_ipos.cast::<f64>() * settings.cell_width();
                                let dist = node_pos - p.x();
                                let node_index = if let Some(rect) = &self.period_bound_rect {
                                    let x_min_index =
                                        (rect.x_min / settings.cell_width()).round() as i64;
                                    let x_max_index =
                                        (rect.x_max / settings.cell_width()).round() as i64;
                                    let y_min_index =
                                        (rect.y_min / settings.cell_width()).round() as i64;
                                    let y_max_index =
                                        (rect.y_max / settings.cell_width()).round() as i64;

                                    let origin = vector![x_min_index, y_min_index];
                                    let node_ipos = node_ipos - origin;
                                    let node_ipos = vector![
                                        (node_ipos.x.rem_euclid(x_max_index - x_min_index)
                                            + origin.x)
                                            as U,
                                        (node_ipos.y.rem_euclid(y_max_index - y_min_index)
                                            + origin.y)
                                            as U
                                    ];

                                    (node_ipos.x, node_ipos.y)
                                } else {
                                    (node_ipos.x as U, node_ipos.y as U)
                                };

                                let params = {
                                    let index = node_index;
                                    if !nodes.contains_key(&index) {
                                        let params = LsmpsParams {
                                            m: Matrix3::<f64>::zeros(),
                                            f_vel: Matrix3x2::<f64>::zeros(),
                                            f_stress: Matrix3::<f64>::zeros(),
                                        };
                                        nodes.insert(index, params);
                                    }

                                    nodes.get_mut(&node_index).unwrap()
                                };
                                let r_ij = dist / rs;
                                let poly_r_ij = poly(r_ij);
                                let weight = weight_function()(
                                    dist.x / settings.cell_width(),
                                    dist.y / settings.cell_width(),
                                );

                                params.m += weight * poly_r_ij * poly_r_ij.transpose();
                                params.f_vel += weight * poly_r_ij.kronecker(&p.v.transpose());
                                let stress =
                                    vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                                params.f_stress +=
                                    weight * poly_r_ij.kronecker(&stress.transpose());
                            });
                    }
                }

                self.grid.par_iter_mut().for_each(|node| {
                    if !nodes.contains_key(&node.index) {
                        return;
                    }
                    let params = nodes.get(&node.index).unwrap();
                    let m_inverse = params.m.try_inverse().unwrap();

                    {
                        let res = scale * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                    }

                    {
                        let res = scale * m_inverse * params.f_stress;
                        node.force[0] = res[(1, 0)] + res[(2, 1)];
                        node.force[1] = res[(1, 1)] + res[(2, 2)];
                    }
                });
            }
            P2GSchemeType::LsmpsOnlyForce => {
                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
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

                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let mass_contrib = node.weight * p.mass;
                        node.node.mass += mass_contrib;
                    }
                }

                let mut nodes = HashMap::new();

                for p in self.particles.iter_mut() {
                    let stress = {
                        let (density, volume) = calc_density_and_volume(
                            settings,
                            p,
                            &self.grid,
                            &self.period_bounds,
                            &self.period_bound_rect,
                        );

                        let mut pressure = 0.;
                        if settings.c != 0. && settings.eos_power != 0. {
                            pressure = settings.rho_0 * settings.c * settings.c
                                / settings.eos_power
                                * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                            if pressure < 0. {
                                pressure = 0.;
                            }
                        }

                        p.pressure = pressure;
                        let pressure = pressure;

                        let dudv = p.c;
                        let strain = dudv;
                        let viscosity_term =
                            settings.dynamic_viscosity * (strain + strain.transpose());

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
                                    p.x().map(|x| (x / settings.cell_width()).floor() as i64)
                                        + vector![gx, gy];
                                let node_pos = node_ipos.cast::<f64>() * settings.cell_width();
                                let dist = node_pos - p.x();
                                if dist.norm() > re {
                                    return;
                                }
                                let node_index = if let Some(rect) = &self.period_bound_rect {
                                    let x_min_index =
                                        (rect.x_min / settings.cell_width()).round() as i64;
                                    let x_max_index =
                                        (rect.x_max / settings.cell_width()).round() as i64;
                                    let y_min_index =
                                        (rect.y_min / settings.cell_width()).round() as i64;
                                    let y_max_index =
                                        (rect.y_max / settings.cell_width()).round() as i64;

                                    let origin = vector![x_min_index, y_min_index];
                                    let node_ipos = node_ipos - origin;
                                    let node_ipos = vector![
                                        (node_ipos.x.rem_euclid(x_max_index - x_min_index)
                                            + origin.x)
                                            as U,
                                        (node_ipos.y.rem_euclid(y_max_index - y_min_index)
                                            + origin.y)
                                            as U
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
                                let r_ij = dist / rs;
                                let poly_r_ij = poly(r_ij);
                                let weight = weight_function()(
                                    dist.x / settings.cell_width(),
                                    dist.y / settings.cell_width(),
                                );
                                // let weight = (1. - (dist / re).norm()).powi(2);

                                params.m += weight * poly_r_ij * poly_r_ij.transpose();
                                let stress =
                                    vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                                params.f_stress +=
                                    weight * poly_r_ij.kronecker(&stress.transpose());
                            });
                    }
                }

                self.grid.par_iter_mut().for_each(|node| {
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
            P2GSchemeType::CompactLsmps => {
                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let mass_contrib = node.weight * p.mass;
                        node.node.mass += mass_contrib;
                    }
                }

                fn factorial(num: usize) -> f64 {
                    match num {
                        0 | 1 => 1.,
                        _ => factorial(num - 1) * num as f64,
                    }
                }

                fn C(bx: usize, by: usize, p: usize, q: usize) -> f64 {
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

                fn S(ax: usize, ay: usize, rs: f64, p: usize, q: usize) -> f64 {
                    let a = ax + ay;
                    let sum = [(0, 0), (1, 0), (0, 1)]
                        .map(|(bx, by)| {
                            let b = bx + by;
                            if b > q || bx > ax || by > ay {
                                0.
                            } else {
                                C(bx, by, p, q) / (factorial(ax - bx) * factorial(ay - by))
                            }
                        })
                        .iter()
                        .sum::<f64>();
                    1. / (sum * rs.powi(a as i32))
                }

                fn poly(r: Vector2<f64>) -> Vector6<f64> {
                    vector![1., r.x, r.y, r.x * r.x, r.x * r.y, r.y * r.y]
                }

                let re = settings.cell_width() * 3.;
                let rs = settings.cell_width();
                let scale_stress = Matrix6::<f64>::from_diagonal(&vector![
                    S(0, 0, rs, 2, 0),
                    S(1, 0, rs, 2, 0),
                    S(0, 1, rs, 2, 0),
                    S(2, 0, rs, 2, 0),
                    S(1, 1, rs, 2, 0),
                    S(0, 2, rs, 2, 0)
                ]);
                let scale_vel = Matrix6::<f64>::from_diagonal(&vector![
                    S(0, 0, rs, 2, 1),
                    S(1, 0, rs, 2, 1),
                    S(0, 1, rs, 2, 1),
                    S(2, 0, rs, 2, 1),
                    S(1, 1, rs, 2, 1),
                    S(0, 2, rs, 2, 1)
                ]);

                struct LsmpsParams {
                    m: Matrix6<f64>,
                    f_vel: Matrix6x2<f64>,
                    f_stress: Matrix6x3<f64>,
                }

                let mut nodes = HashMap::new();

                for p in self.particles.iter_mut() {
                    let stress = {
                        let (density, volume) = calc_density_and_volume(
                            settings,
                            p,
                            &self.grid,
                            &self.period_bounds,
                            &self.period_bound_rect,
                        );

                        let mut pressure = 0.;
                        if settings.c != 0. && settings.eos_power != 0. {
                            pressure = settings.rho_0 * settings.c * settings.c
                                / settings.eos_power
                                * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                            if pressure < 0. {
                                pressure = 0.;
                            }
                        }

                        // pressure = {
                        //     let PI = std::f64::consts::PI;
                        //     let L = 1.;
                        //     let rho = 1.;
                        //     let U = 1.;
                        //     let nu = 1e-2;

                        //     let (x, y) = (p.x.x - 5., p.x.y - 5.);

                        //     rho * U * U / 4.
                        //         * f64::exp(
                        //             -4. * PI * PI * (self.steps as f64) * settings.dt * nu
                        //                 / (L * L),
                        //         )
                        //         * (f64::cos(2. * PI * x / L) + f64::cos(2. * PI * y / L))
                        // };

                        p.pressure = pressure;
                        let pressure = pressure;

                        let dudv = p.c;
                        let strain = dudv;
                        let viscosity_term =
                            settings.dynamic_viscosity * (strain + strain.transpose());

                        (-pressure * Matrix2f::identity() + viscosity_term)
                    };

                    for node in NodeIterator::new(
                        settings,
                        &self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let params = {
                            let index = node.node.index;
                            if !nodes.contains_key(&index) {
                                let params = LsmpsParams {
                                    m: Matrix6::<f64>::zeros(),
                                    f_vel: Matrix6x2::<f64>::zeros(),
                                    f_stress: Matrix6x3::<f64>::zeros(),
                                };
                                nodes.insert(index, params);
                            }

                            nodes.get_mut(&node.node.index).unwrap()
                        };

                        let r_ij = -node.dist / rs;
                        let poly_r_ij = poly(r_ij);
                        let weight = node.weight;

                        params.m += weight * poly_r_ij * poly_r_ij.transpose();

                        params.f_vel +=
                            weight * poly_r_ij.kronecker(&p.v.transpose()) * C(0, 0, 2, 1) as f64;
                        params.f_vel += weight
                            * poly_r_ij.kronecker(&p.c.column(0).transpose())
                            * C(1, 0, 2, 1) as f64
                            * -node.dist.x;
                        params.f_vel += weight
                            * poly_r_ij.kronecker(&p.c.column(1).transpose())
                            * C(0, 1, 2, 1) as f64
                            * -node.dist.y;

                        let stress = vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                        params.f_stress += weight
                            * poly_r_ij.kronecker(&stress.transpose())
                            * C(0, 0, 2, 0) as f64;
                    }
                }

                self.grid.par_iter_mut().for_each(|node| {
                    if !nodes.contains_key(&node.index) {
                        return;
                    }
                    let params = nodes.get(&node.index).unwrap();
                    let m_inverse = params.m.try_inverse().unwrap();

                    {
                        let res = scale_vel * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                    }

                    {
                        let res = scale_stress * m_inverse * params.f_stress;
                        node.force[0] = res[(1, 0)] + res[(2, 1)];
                        node.force[1] = res[(1, 1)] + res[(2, 2)];
                    }
                });
            }
            P2GSchemeType::CompactLsmpsLinear => {
                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let mass_contrib = node.weight * p.mass;
                        node.node.mass += mass_contrib;
                    }
                }

                fn factorial(num: usize) -> f64 {
                    match num {
                        0 | 1 => 1.,
                        _ => factorial(num - 1) * num as f64,
                    }
                }

                fn C(bx: usize, by: usize, p: usize, q: usize) -> f64 {
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

                fn S(ax: usize, ay: usize, rs: f64, p: usize, q: usize) -> f64 {
                    let a = ax + ay;
                    let sum = [(0, 0), (1, 0), (0, 1)]
                        .map(|(bx, by)| {
                            let b = bx + by;
                            if b > q || bx > ax || by > ay {
                                0.
                            } else {
                                C(bx, by, p, q) / (factorial(ax - bx) * factorial(ay - by))
                            }
                        })
                        .iter()
                        .sum::<f64>();
                    1. / (sum * rs.powi(a as i32))
                }

                fn poly(r: Vector2<f64>) -> Vector3<f64> {
                    vector![1., r.x, r.y]
                }

                let re = settings.cell_width() * 3.;
                let rs = settings.cell_width();
                let scale_stress = Matrix3::<f64>::from_diagonal(&vector![
                    S(0, 0, rs, 1, 0),
                    S(1, 0, rs, 1, 0),
                    S(0, 1, rs, 1, 0)
                ]);
                let scale_vel = Matrix3::<f64>::from_diagonal(&vector![
                    S(0, 0, rs, 1, 1),
                    S(1, 0, rs, 1, 1),
                    S(0, 1, rs, 1, 1)
                ]);

                struct LsmpsParams {
                    m: Matrix3<f64>,
                    f_vel: Matrix3x2<f64>,
                    f_stress: Matrix3<f64>,
                }

                let mut nodes = HashMap::new();

                for p in self.particles.iter_mut() {
                    let stress = {
                        let (density, volume) = calc_density_and_volume(
                            settings,
                            p,
                            &self.grid,
                            &self.period_bounds,
                            &self.period_bound_rect,
                        );

                        let mut pressure = 0.;
                        if settings.c != 0. && settings.eos_power != 0. {
                            pressure = settings.rho_0 * settings.c * settings.c
                                / settings.eos_power
                                * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                            if pressure < 0. {
                                pressure = 0.;
                            }
                        }

                        // pressure = {
                        //     let PI = std::f64::consts::PI;
                        //     let L = 1.;
                        //     let rho = 1.;
                        //     let U = 1.;
                        //     let nu = 1e-2;

                        //     let (x, y) = (p.x.x - 5., p.x.y - 5.);

                        //     rho * U * U / 4.
                        //         * f64::exp(
                        //             -4. * PI * PI * (self.steps as f64) * settings.dt * nu
                        //                 / (L * L),
                        //         )
                        //         * (f64::cos(2. * PI * x / L) + f64::cos(2. * PI * y / L))
                        // };

                        p.pressure = pressure;
                        let pressure = pressure;

                        let dudv = p.c;
                        let strain = dudv;
                        let viscosity_term =
                            settings.dynamic_viscosity * (strain + strain.transpose());

                        (-pressure * Matrix2f::identity() + viscosity_term)
                    };

                    for node in NodeIterator::new(
                        settings,
                        &self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let params = {
                            let index = node.node.index;
                            if !nodes.contains_key(&index) {
                                let params = LsmpsParams {
                                    m: Matrix3::<f64>::zeros(),
                                    f_vel: Matrix3x2::<f64>::zeros(),
                                    f_stress: Matrix3::<f64>::zeros(),
                                };
                                nodes.insert(index, params);
                            }

                            nodes.get_mut(&node.node.index).unwrap()
                        };

                        let r_ij = -node.dist / rs;
                        let poly_r_ij = poly(r_ij);
                        let weight = node.weight;

                        params.m += weight * poly_r_ij * poly_r_ij.transpose();

                        params.f_vel +=
                            weight * poly_r_ij.kronecker(&p.v.transpose()) * C(0, 0, 1, 1);
                        params.f_vel += weight
                            * poly_r_ij.kronecker(&p.c.column(0).transpose())
                            * C(1, 0, 1, 1)
                            * -node.dist.x;
                        params.f_vel += weight
                            * poly_r_ij.kronecker(&p.c.column(1).transpose())
                            * C(0, 1, 1, 1)
                            * -node.dist.y;

                        let stress = vector![stress[(0, 0)], stress[(0, 1)], stress[(1, 1)]];
                        params.f_stress +=
                            weight * poly_r_ij.kronecker(&stress.transpose()) * C(0, 0, 1, 0);
                    }
                }

                self.grid.par_iter_mut().for_each(|node| {
                    if !nodes.contains_key(&node.index) {
                        return;
                    }
                    let params = nodes.get(&node.index).unwrap();
                    let m_inverse = params.m.try_inverse().unwrap();

                    {
                        let res = scale_vel * m_inverse * params.f_vel;
                        node.v = res.row(0).transpose();
                    }

                    {
                        let res = scale_stress * m_inverse * params.f_stress;
                        node.force[0] = res[(1, 0)] + res[(2, 1)];
                        node.force[1] = res[(1, 1)] + res[(2, 2)];
                    }
                });
            }
            P2GSchemeType::CompactOnlyVelocity => {
                for p in self.particles.iter() {
                    for node in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let mass_contrib = node.weight * p.mass;
                        node.node.mass += mass_contrib;
                    }
                }

                for p in self.particles.iter_mut() {
                    let (density, volume) = calc_density_and_volume(
                        settings,
                        p,
                        &self.grid,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    );

                    let mut pressure = 0.;
                    if settings.c != 0. && settings.eos_power != 0. {
                        pressure = settings.rho_0 * settings.c * settings.c / settings.eos_power
                            * ((density / settings.rho_0).powf(settings.eos_power) - 1.);
                        if pressure < 0. {
                            pressure = 0.;
                        }
                    }

                    p.pressure = pressure;
                    let pressure = pressure;

                    let dudv = p.c;
                    let strain = dudv;
                    let viscosity_term = settings.dynamic_viscosity * (strain + strain.transpose());
                    let stress = -pressure * Matrix2f::identity() + viscosity_term;
                    let eq_16_term_0 =
                        -volume * 4. / (settings.cell_width() * settings.cell_width()) * stress;

                    for n in NodeMutIterator::new(
                        settings,
                        &mut self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        n.node.force += eq_16_term_0 * n.weight * n.dist;
                    }
                }

                fn factorial(num: usize) -> f64 {
                    match num {
                        0 | 1 => 1.,
                        _ => factorial(num - 1) * num as f64,
                    }
                }

                fn C(bx: usize, by: usize, p: usize, q: usize) -> f64 {
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

                fn S(ax: usize, ay: usize, rs: f64, p: usize, q: usize) -> f64 {
                    let a = ax + ay;
                    let sum = [(0, 0), (1, 0), (0, 1)]
                        .map(|(bx, by)| {
                            let b = bx + by;
                            if b > q || bx > ax || by > ay {
                                0.
                            } else {
                                C(bx, by, p, q) / (factorial(ax - bx) * factorial(ay - by))
                            }
                        })
                        .iter()
                        .sum::<f64>();
                    1. / (sum * rs.powi(a as i32))
                }

                fn poly(r: Vector2<f64>) -> Vector6<f64> {
                    vector![1., r.x, r.y, r.x * r.x, r.x * r.y, r.y * r.y]
                }

                let re = settings.cell_width() * 3.;
                let rs = settings.cell_width();
                let scale_vel = Matrix6::<f64>::from_diagonal(&vector![
                    S(0, 0, rs, 2, 1),
                    S(1, 0, rs, 2, 1),
                    S(0, 1, rs, 2, 1),
                    S(2, 0, rs, 2, 1),
                    S(1, 1, rs, 2, 1),
                    S(0, 2, rs, 2, 1)
                ]);

                struct LsmpsParams {
                    m: Matrix6<f64>,
                    f_vel: Matrix6x2<f64>,
                }

                let mut nodes = HashMap::new();

                for p in self.particles.iter_mut() {
                    for node in NodeIterator::new(
                        settings,
                        &self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let params = {
                            let index = node.node.index;
                            if !nodes.contains_key(&index) {
                                let params = LsmpsParams {
                                    m: Matrix6::<f64>::zeros(),
                                    f_vel: Matrix6x2::<f64>::zeros(),
                                };
                                nodes.insert(index, params);
                            }

                            nodes.get_mut(&node.node.index).unwrap()
                        };

                        let r_ij = node.dist / rs;
                        let poly_r_ij = poly(r_ij);
                        let weight = node.weight;

                        params.m += weight * poly_r_ij * poly_r_ij.transpose();

                        params.f_vel +=
                            weight * poly_r_ij.kronecker(&p.v.transpose()) * C(0, 0, 2, 1) as f64;
                        params.f_vel += weight
                            * poly_r_ij.kronecker(&p.c.column(0).transpose())
                            * C(1, 0, 2, 1) as f64
                            * -node.dist.x;
                        params.f_vel += weight
                            * poly_r_ij.kronecker(&p.c.column(1).transpose())
                            * C(0, 1, 2, 1) as f64
                            * -node.dist.y;
                    }
                }

                self.grid.par_iter_mut().for_each(|node| {
                    if !nodes.contains_key(&node.index) {
                        return;
                    }
                    let params = nodes.get(&node.index).unwrap();
                    let m_inverse = params.m.try_inverse().unwrap();

                    let res = scale_vel * m_inverse * params.f_vel;
                    node.v = res.row(0).transpose();
                });
            }
        }
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
                | P2GSchemeType::CompactLsmps
                | P2GSchemeType::CompactLsmpsLinear => {
                    n.v_star = n.v + settings.dt * (vector![0., settings.gravity] + n.force);
                }
                P2GSchemeType::CompactOnlyVelocity => {
                    n.v_star =
                        n.v + settings.dt * (vector![0., settings.gravity] + n.force / n.mass);
                }
                P2GSchemeType::LsmpsOnlyForce => {
                    n.v /= n.mass;
                    n.v_star = n.v + settings.dt * (vector![0., settings.gravity] + n.force);
                }
                _ => {
                    n.v /= n.mass;
                    n.v_star =
                        n.v + settings.dt * (vector![0., settings.gravity] + n.force / n.mass);
                }
            };

            let node_pos = calc_node_pos(settings, i);
            for b in self.slip_bounds.iter() {
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
    }

    pub fn g2p(&mut self, settings: &Settings) {
        self.particles.par_iter_mut().for_each(|p| {
            match settings.g2p_scheme {
                G2PSchemeType::MLSMPM => {
                    let p_v_t = p.v;
                    p.v = Vector2f::zeros();
                    p.c = Matrix2f::zeros();

                    for n in NodeIterator::new(
                        settings,
                        &self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
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
                }
                G2PSchemeType::LSMPS => {
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
                        &self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let r_ij = n.dist / rs;
                        let poly_r_ij = poly(r_ij);
                        let weight = n.weight;

                        params.m += weight * poly_r_ij * poly_r_ij.transpose();
                        params.f_vel += weight * poly_r_ij.kronecker(&n.node.v_star.transpose());
                    }

                    if let Ok(m_inverted) = params.m.pseudo_inverse(1e-15) {
                        let res = scale * m_inverted * params.f_vel;
                        p.v = res.row(0).transpose();
                        p.x += res.row(0).transpose() * settings.dt;
                        p.c = Matrix2::new(res.m21, res.m31, res.m22, res.m32);
                    }
                }
                G2PSchemeType::LsmpsLinear => {
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
                        &self.grid,
                        p,
                        &self.period_bounds,
                        &self.period_bound_rect,
                    ) {
                        let r_ij = n.dist / rs;
                        let poly_r_ij = poly(r_ij);
                        let weight = n.weight;

                        params.m += weight * poly_r_ij * poly_r_ij.transpose();
                        params.f_vel += weight * poly_r_ij.kronecker(&n.node.v_star.transpose());
                    }

                    if let Ok(m_inverted) = params.m.pseudo_inverse(1e-15) {
                        let res = scale * m_inverted * params.f_vel;
                        //println!("{}", m_inverted);
                        p.v = res.row(0).transpose();
                        p.x += res.row(0).transpose() * settings.dt;
                        p.c = Matrix2::new(res.m21, res.m31, res.m22, res.m32);
                    }
                }
            }

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
}
