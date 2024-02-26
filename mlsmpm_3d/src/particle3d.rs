use crate::*;

#[derive(Debug, Clone)]
pub struct Particle {
    pub x: Vector3f,
    pub v: Vector3f,
    pub c: Matrix3f,
    pub mass: f64,
    pub pressure: f64,
}

impl Particle {
    pub fn new(pos: Vector3f) -> Particle {
        Particle {
            x: pos,
            v: Vector3::zeros(),
            c: Matrix3f::zeros(),
            mass: 0.,
            pressure: 0.,
        }
    }

    pub fn x(&self) -> &Vector3f {
        &self.x
    }

    pub fn v(&self) -> &Vector3f {
        &self.v
    }

    pub fn pressure(&self) -> f64 {
        self.pressure
    }

    pub fn v_norm(&self) -> f64 {
        self.v.norm()
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
