use crate::*;

#[derive(Debug, Clone)]
pub struct Node {
    pub(super) v: Vector3f,
    pub(super) v_star: Vector3f,
    pub(super) force: Vector3f,
    pub(super) mass: f64,
    pub(super) basis_function: fn(f64) -> f64,
    pub(super) node_type: NodeType,
}

#[derive(Debug, PartialEq, Clone)]
pub enum NodeType {
    Normal,
    LeftWall,
    LeftHalf,
    RightWall,
    RightHalf,
}

impl Node {
    pub fn new(node_type: NodeType) -> Node {
        let basis = match node_type {
            NodeType::Normal => quadratic_b_spline,
            NodeType::LeftWall => quadratic_b_spline,
            NodeType::LeftHalf => quadratic_b_spline,
            NodeType::RightWall => quadratic_b_spline,
            NodeType::RightHalf => quadratic_b_spline,
        };

        Node {
            v: Vector3::zeros(),
            v_star: Vector3::zeros(),
            force: Vector3::zeros(),
            mass: 0.0,
            basis_function: basis,
            node_type,
        }
    }

    pub fn reset(&mut self) {
        self.v = Vector3::zeros();
        self.v_star = Vector3::zeros();
        self.force = Vector3::zeros();
        self.mass = 0.0;
    }

    pub fn calc_weight(&self, x: Vector3f, h: f64) -> f64 {
        let wx = (self.basis_function)(x.x / h);
        let wy = (self.basis_function)(x.y / h);
        let wz = (self.basis_function)(x.z / h);

        wx * wy * wz
    }

    pub fn formatted_list(&self) -> [String; 7] {
        [
            self.v.x.to_string(),
            self.v.y.to_string(),
            self.v_star.x.to_string(),
            self.v_star.y.to_string(),
            self.force.x.to_string(),
            self.force.y.to_string(),
            self.mass.to_string(),
        ]
    }
}

fn quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1.5, -0.5, 0.5, 1.5], x)
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

fn left_half_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1., -0.5, 0.5, 1.5], x)
}

fn left_wall_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![0., 0., 0.5, 1.5], x)
        + b_spline_basis(0, 2, &vec![0., 0., 0., 0.5], x)
}

fn right_half_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1.5, -0.5, 0.5, 1.], x)
}

fn right_wall_quadratic_kernel(x: f64) -> f64 {
    b_spline_basis(0, 2, &vec![-1.5, -0.5, 0., 0.], x)
        + b_spline_basis(0, 2, &vec![-0.5, 0., 0., 0.], x)
}

fn b_spline_basis(j: usize, k: usize, knot: &Vec<f64>, t: f64) -> f64 {
    if knot.len() as i32 - k as i32 - 2 < 0 || j > knot.len() - 2 {
        panic!();
    }

    if k == 0 {
        if knot[j] <= t && t < knot[j + 1] {
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

        w(j, k, knot, t) * b_spline_basis(j, k - 1, knot, t)
            + (1. - w(j + 1, k, knot, t)) * b_spline_basis(j + 1, k - 1, knot, t)
    }
}
