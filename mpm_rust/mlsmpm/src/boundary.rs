use crate::*;

#[derive(Debug)]
pub enum Direction {
    X,
    Y,
}

#[derive(Debug)]
pub struct BoundaryLine<T: std::ops::Sub<Output = T>> {
    pub(super) value: T,
    pub(super) lower: bool,
}

impl<T: std::ops::Sub<Output = T> + std::ops::Add<Output = T> + std::marker::Copy> BoundaryLine<T> {
    pub fn new(value: T, lower: bool) -> Self {
        Self { value, lower }
    }

    pub fn calc_excess(&self, target: T) -> T {
        if self.lower {
            self.value - target
        } else {
            target - self.value
        }
    }

    pub fn plus_excess(&self, excess: T) -> T {
        if self.lower {
            self.value + excess
        } else {
            self.value - excess
        }
    }
}

#[derive(Debug)]
pub struct SlipBoundary {
    pub(super) line: BoundaryLine<f64>,
    pub(super) direction: Direction,
}

impl SlipBoundary {
    pub fn new(value: f64, dir: Direction, lower: bool) -> SlipBoundary {
        SlipBoundary {
            line: BoundaryLine::<f64> { value, lower },
            direction: dir,
        }
    }
}

#[derive(Debug)]
pub struct PeriodicBoundary {
    pub(super) a: BoundaryLine<f64>,
    pub(super) b: BoundaryLine<f64>,
    pub(super) direction: Direction,
}

impl PeriodicBoundary {
    pub fn new(a: BoundaryLine<f64>, b: BoundaryLine<f64>, dir: Direction) -> PeriodicBoundary {
        PeriodicBoundary {
            a,
            b,
            direction: dir,
        }
    }
}
