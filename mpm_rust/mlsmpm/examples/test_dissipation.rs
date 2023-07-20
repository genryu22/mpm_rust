use chrono::Local;
use nalgebra::*;
use rand::Rng;
use rayon::prelude::{
    IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};
use std::{error::Error, fs, path::Path};

struct Particle {
    x: f64,
    v: f64,
    c: f64,
}

struct Node {
    i: usize,
    v: f64,
}

fn init(
    space_size: usize,
    grid_size: usize,
    init_vel: fn(f64) -> f64,
    init_vel_grad: fn(f64) -> f64,
) -> (Vec<Particle>, Vec<Node>, f64) {
    let grid = (0..grid_size + 1)
        .map(|index| Node { i: index, v: 0. })
        .collect::<Vec<_>>();

    let cell_width = space_size as f64 / grid_size as f64;

    let particles = (0..2 * grid_size + 1)
        .map(|i| {
            let x = i as f64 * cell_width * 0.5;
            Particle {
                x,
                v: init_vel(x),
                c: init_vel_grad(x),
            }
        })
        .collect::<Vec<_>>();

    (particles, grid, cell_width)
}

fn sin_vel(x: f64) -> f64 {
    let PI = std::f64::consts::PI;
    f64::sin(PI * x * 4.)
}

fn sin_vel_grad(x: f64) -> f64 {
    let PI = std::f64::consts::PI;
    4. * PI * f64::cos(PI * x * 4.)
}

fn main() -> Result<(), Box<dyn Error>> {
    let folder = Path::new("exp_dissipation");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let cases: [(
        fn(&Vec<Particle>, &mut Vec<Node>, f64, usize),
        fn(&mut Vec<Particle>, &Vec<Node>, f64, usize),
        &str,
    ); 3] = [
        (mlsmpm_p2g, mlsmpm_g2p, "mlsmpm"),
        (lsmps_p2g, lsmps_g2p, "lsmps"),
        (compact_lsmps_p2g, lsmps_g2p, "compact_lsmps"),
    ];

    let results = cases
        .par_iter()
        .map(|(p2g, g2p, name)| {
            let grid_size = 300;
            let (mut particles, mut grid, cell_width) = init(1, grid_size, sin_vel, sin_vel_grad);
            for _ in 0..100 {
                p2g(&particles, &mut grid, cell_width, grid_size);
                particles.par_iter_mut().for_each(|p| p.v = 0.);
                g2p(&mut particles, &grid, cell_width, grid_size);
                grid.par_iter_mut().for_each(|n| n.v = 0.);
            }
            (name, particles, grid)
        })
        .collect::<Vec<_>>();

    let dt = Local::now();
    for r in results {
        let mut writer = csv::Writer::from_path(folder.join(format!(
            "{:}_{:}.csv",
            r.0,
            dt.format("%Y%m%d_%Hh%Mm%Ss")
        )))?;

        writer.write_record(&["x", "v", "true_v"])?;
        let particles = r.1;
        for p in particles {
            writer.write_record(&[p.x.to_string(), p.v.to_string(), sin_vel(p.x).to_string()])?;
        }
        writer.flush()?;
    }

    Ok(())
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

fn mlsmpm_p2g(particles: &Vec<Particle>, grid: &mut Vec<Node>, cell_width: f64, grid_size: usize) {
    particles.iter().for_each(|p| {
        let e = 3;
        (-e..=e)
            .map(|gx| {
                let base_node = (p.x / cell_width).floor() as i32;
                (
                    (base_node + gx) as f64 * cell_width,
                    (base_node + gx).rem_euclid(grid_size as i32 + 1) as usize,
                )
            })
            .for_each(|(n_pos, n_index)| {
                let n = &mut grid[n_index];
                let dist = n_pos - p.x;
                let weight = quadratic_b_spline(dist / cell_width);
                n.v += 0.5 * weight * p.v;
            });
    });
}

fn mlsmpm_g2p(particles: &mut Vec<Particle>, grid: &Vec<Node>, cell_width: f64, grid_size: usize) {
    particles.iter_mut().for_each(|p| {
        let e = 3;
        (-e..=e)
            .map(|gx| {
                let base_node = (p.x / cell_width).floor() as i32;
                (
                    (base_node + gx) as f64 * cell_width,
                    (base_node + gx).rem_euclid(grid_size as i32 + 1) as usize,
                )
            })
            .for_each(|(n_pos, n_index)| {
                let n = &grid[n_index];
                let dist = n_pos - p.x;
                let weight = quadratic_b_spline(dist / cell_width);
                p.v += weight * n.v;
            });
    });
}

fn lsmps_p2g(particles: &Vec<Particle>, grid: &mut Vec<Node>, cell_width: f64, grid_size: usize) {
    fn poly(r: f64) -> Vector3<f64> {
        vector![1., r, r * r]
    }

    let rs = cell_width;
    let scale = Matrix3::<f64>::from_diagonal(&vector![1., 1. / rs, 2. / rs / rs]);

    struct LsmpsParams {
        m: Matrix3<f64>,
        f_vel: Vector3<f64>,
    }

    let params = particles
        .par_iter()
        .flat_map(|p| {
            let e = 10;
            (-e..=e)
                .map(|gx| {
                    let base_node = (p.x / cell_width).floor() as i32;
                    (
                        (base_node + gx) as f64 * cell_width,
                        (base_node + gx).rem_euclid(grid_size as i32 + 1) as usize,
                    )
                })
                .filter(|(n_pos, _)| (n_pos - p.x).abs() <= e as f64 * cell_width)
                .map(|(n_pos, n_index)| {
                    let dist = n_pos - p.x;
                    let r_ij = -dist / rs;
                    let poly_r_ij = poly(r_ij);
                    //let weight = quadratic_b_spline(dist / cell_width);
                    let weight = (1. - (dist / (e as f64 * cell_width)).abs()).powi(2);

                    (
                        n_index,
                        LsmpsParams {
                            m: weight * poly_r_ij * poly_r_ij.transpose(),
                            f_vel: weight * poly_r_ij.kronecker(&Vector1::new(p.v)),
                        },
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut params_sums = {
        let mut params_sums = vec![];
        for _ in 0..grid.len() {
            params_sums.push(LsmpsParams {
                m: Matrix3::zeros(),
                f_vel: Vector3::zeros(),
            });
        }
        params_sums
    };
    for (n_index, params) in params {
        params_sums[n_index].m += params.m;
        params_sums[n_index].f_vel += params.f_vel;
    }

    for (index, vel) in params_sums
        .par_iter_mut()
        .map(|params| {
            if let Some(m_inverse) = params.m.try_inverse() {
                let res = scale * m_inverse * params.f_vel;
                res.row(0).transpose().x
            } else {
                0.
            }
        })
        .enumerate()
        .collect::<Vec<_>>()
    {
        grid[index].v = vel;
    }
}

fn lsmps_g2p(particles: &mut Vec<Particle>, grid: &Vec<Node>, cell_width: f64, grid_size: usize) {
    fn poly(r: f64) -> Vector3<f64> {
        vector![1., r, r * r]
    }

    let rs = cell_width;
    let scale = Matrix3::<f64>::from_diagonal(&vector![1., 1. / rs, 2. / rs / rs]);

    struct LsmpsParams {
        m: Matrix3<f64>,
        f_vel: Vector3<f64>,
    }

    impl std::iter::Sum for LsmpsParams {
        fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
            iter.reduce(|mut acc, e| {
                acc.m += e.m;
                acc.f_vel += e.f_vel;
                acc
            })
            .unwrap_or(LsmpsParams {
                m: Matrix3::zeros(),
                f_vel: Vector3::zeros(),
            })
        }
    }

    particles.iter_mut().for_each(|p| {
        let e = 10;
        let params = (-e..=e)
            .map(|gx| {
                let base_node = (p.x / cell_width).floor() as i32;
                (
                    (base_node + gx) as f64 * cell_width,
                    (base_node + gx).rem_euclid(grid_size as i32 + 1) as usize,
                )
            })
            .filter(|(n_pos, _)| (n_pos - p.x).abs() <= e as f64 * cell_width)
            .map(|(n_pos, n_index)| {
                let dist = n_pos - p.x;
                let r_ij = dist / rs;
                let poly_r_ij = poly(r_ij);
                //let weight = quadratic_b_spline(dist / cell_width);
                let weight = (1. - (dist / (e as f64 * cell_width)).abs()).powi(2);

                LsmpsParams {
                    m: weight * poly_r_ij * poly_r_ij.transpose(),
                    f_vel: weight * poly_r_ij.kronecker(&Vector1::new(grid[n_index].v)),
                }
            })
            .sum::<LsmpsParams>();

        if let Some(m_inverted) = params.m.try_inverse() {
            let res = scale * m_inverted * params.f_vel;
            p.v = res.row(0).x;
            p.c = res.row(1).x;
        }
    });
}

fn compact_lsmps_p2g(
    particles: &Vec<Particle>,
    grid: &mut Vec<Node>,
    cell_width: f64,
    grid_size: usize,
) {
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

    fn poly(r: f64) -> Vector3<f64> {
        vector![1., r, r * r]
    }

    let rs = cell_width;
    let scale = Matrix3::<f64>::from_diagonal(&vector![
        S(0, 0, rs, 2, 1),
        S(1, 0, rs, 2, 1),
        S(2, 0, rs, 2, 1)
    ]);

    struct LsmpsParams {
        m: Matrix3<f64>,
        f_vel: Vector3<f64>,
    }

    let params = particles
        .par_iter()
        .flat_map(|p| {
            let e = 10;
            (-e..=e)
                .map(|gx| {
                    let base_node = (p.x / cell_width).floor() as i32;
                    (
                        (base_node + gx) as f64 * cell_width,
                        (base_node + gx).rem_euclid(grid_size as i32 + 1) as usize,
                    )
                })
                .filter(|(n_pos, _)| (n_pos - p.x).abs() <= e as f64 * cell_width)
                .map(|(n_pos, n_index)| {
                    let dist = n_pos - p.x;
                    let r_ij = -dist / rs;
                    let poly_r_ij = poly(r_ij);
                    //let weight = quadratic_b_spline(dist / cell_width);
                    let weight = (1. - (dist / (e as f64 * cell_width)).abs()).powi(2);

                    (
                        n_index,
                        LsmpsParams {
                            m: weight * poly_r_ij * poly_r_ij.transpose(),
                            f_vel: weight * poly_r_ij.kronecker(&Vector1::new(p.v))
                                + weight
                                    * poly_r_ij.kronecker(&Vector1::new(p.c))
                                    * C(1, 0, 2, 1) as f64
                                    * -dist,
                        },
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut params_sums = {
        let mut params_sums = vec![];
        for _ in 0..grid.len() {
            params_sums.push(LsmpsParams {
                m: Matrix3::zeros(),
                f_vel: Vector3::zeros(),
            });
        }
        params_sums
    };
    for (n_index, params) in params {
        params_sums[n_index].m += params.m;
        params_sums[n_index].f_vel += params.f_vel;
    }

    for (index, vel) in params_sums
        .par_iter_mut()
        .map(|params| {
            if let Some(m_inverse) = params.m.try_inverse() {
                let res = scale * m_inverse * params.f_vel;
                res.row(0).transpose().x
            } else {
                0.
            }
        })
        .enumerate()
        .collect::<Vec<_>>()
    {
        grid[index].v = vel;
    }
}
