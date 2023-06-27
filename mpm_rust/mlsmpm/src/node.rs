use crate::*;

#[derive(Debug, Clone)]
pub struct Node {
    pub(super) v: Vector2f,
    pub(super) v_star: Vector2f,
    pub(super) force: Vector2f,
    pub(super) mass: f64,
    pub(super) index: (usize, usize),
}

impl Node {
    pub fn new(index: (usize, usize)) -> Node {
        Node {
            v: Vector2::zeros(),
            v_star: Vector2::zeros(),
            force: Vector2::zeros(),
            mass: 0.0,
            index,
        }
    }

    pub fn reset(&mut self) {
        self.v = Vector2::zeros();
        self.v_star = Vector2::zeros();
        self.force = Vector2::zeros();
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
