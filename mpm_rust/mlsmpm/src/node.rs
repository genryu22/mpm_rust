use crate::*;

#[derive(Debug)]
pub struct Node {
    pub(super) v: Vector2f,
    pub(super) v_star: Vector2f,
    pub(super) force: Vector2f,
    pub(super) mass: f64,
}

impl Node {
    pub fn new() -> Node {
        Node {
            v: Vector2::zeros(),
            v_star: Vector2::zeros(),
            force: Vector2::zeros(),
            mass: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.v = Vector2::zeros();
        self.v_star = Vector2::zeros();
        self.force = Vector2::zeros();
        self.mass = 0.0;
    }
}
