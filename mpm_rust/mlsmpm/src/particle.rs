use crate::*;

#[derive(Debug, Clone)]
pub struct Particle {
    pub(super) x: Vector2f,
    pub(super) v: Vector2f,
    pub(super) c: Matrix2f,
    pub(super) mass: f64,
}

impl Particle {
    pub fn new(pos: Vector2f) -> Particle {
        Particle {
            x: pos,
            v: Vector2::zeros(),
            c: Matrix2f::zeros(),
            mass: 0.0,
        }
    }

    pub fn x(&self) -> &Vector2f {
        &self.x
    }

    pub fn v(&self) -> &Vector2f {
        &self.v
    }
}
