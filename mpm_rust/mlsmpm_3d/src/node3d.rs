use crate::*;

#[derive(Debug)]
pub struct Node {
    pub(super) v: Vector3f,
    pub(super) v_star: Vector3f,
    pub(super) force: Vector3f,
    pub(super) mass: f64,
}

impl Node {
    pub fn new() -> Node {
        Node {
            v: Vector3::zeros(),
            v_star: Vector3::zeros(),
            force: Vector3::zeros(),
            mass: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.v = Vector3::zeros();
        self.v_star = Vector3::zeros();
        self.force = Vector3::zeros();
        self.mass = 0.0;
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
