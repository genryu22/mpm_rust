use crate::*;

#[derive(Debug, Clone)]
pub struct Particle {
    pub(super) x: Vector2f,
    pub(super) v: Vector2f,
    pub(super) c: Matrix2f,
    pub(super) mass: f64,
    pub(super) pressure: f64,
}

impl Particle {
    pub fn new(pos: Vector2f) -> Particle {
        Particle {
            x: pos,
            v: Vector2::zeros(),
            c: Matrix2f::zeros(),
            mass: 0.,
            pressure: 0.,
        }
    }

    pub fn x(&self) -> &Vector2f {
        &self.x
    }

    pub fn v(&self) -> &Vector2f {
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
