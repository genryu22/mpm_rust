use na::Matrix2xX;

use crate::*;

#[derive(Debug, Clone)]
pub struct Particle {
    pub(super) x: Vector2f,
    pub(super) init_x: Vector2f,
    pub(super) v: Vector2f,
    pub(super) c: Matrix2f,
    pub(super) mass: f64,
    pub(super) pressure: f64,
    pub(super) v_lsmps: Matrix2xX<f64>,
}

impl Particle {
    pub fn new(pos: Vector2f) -> Particle {
        Particle {
            x: pos,
            init_x: pos,
            v: Vector2::zeros(),
            c: Matrix2f::zeros(),
            mass: 0.,
            pressure: 0.,
            v_lsmps: Matrix2xX::zeros(0),
        }
    }

    pub fn new_with_mass(pos: Vector2f, mass: f64) -> Particle {
        Particle {
            x: pos,
            init_x: pos,
            v: Vector2::zeros(),
            c: Matrix2f::zeros(),
            mass,
            pressure: 0.,
            v_lsmps: Matrix2xX::zeros(0),
        }
    }

    pub fn new_with_mass_velocity(pos: Vector2f, mass: f64, velocity: Vector2f) -> Particle {
        Particle {
            x: pos,
            init_x: pos,
            v: velocity,
            c: Matrix2f::zeros(),
            mass,
            pressure: 0.,
            v_lsmps: Matrix2xX::zeros(0),
        }
    }

    pub fn new_with_mass_velocity_c(
        pos: Vector2f,
        mass: f64,
        velocity: Vector2f,
        c: Matrix2f,
    ) -> Particle {
        Particle {
            x: pos,
            init_x: pos,
            v: velocity,
            c,
            mass,
            pressure: 0.,
            v_lsmps: Matrix2xX::zeros(0),
        }
    }

    pub fn new_with_mass_velocity_c_lsmps(
        pos: Vector2f,
        mass: f64,
        velocity: Vector2f,
        c: Matrix2f,
        x_lsmps: Matrix2xX<f64>,
    ) -> Particle {
        Particle {
            x: pos,
            init_x: pos,
            v: velocity,
            c,
            mass,
            pressure: 0.,
            v_lsmps: x_lsmps,
        }
    }

    pub fn mass(&self) -> f64 {
        self.mass
    }

    pub fn x(&self) -> &Vector2f {
        &self.x
    }

    pub fn v(&self) -> &Vector2f {
        &self.v
    }

    pub fn c(&self) -> &Matrix2f {
        &self.c
    }

    pub fn pressure(&self) -> f64 {
        self.pressure
    }

    pub fn v_norm(&self) -> f64 {
        self.v.norm()
    }

    pub fn x_lsmps(&self) -> Matrix2xX<f64> {
        self.v_lsmps.clone()
    }

    pub fn formatted_list(&self) -> [String; 6] {
        [
            self.x.x.to_string(),
            self.x.y.to_string(),
            self.v.x.to_string(),
            self.v.y.to_string(),
            self.mass.to_string(),
            self.pressure.to_string(),
        ]
    }
}
