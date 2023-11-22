

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
    pub(super) fixed: bool,
    pub(super) no_slip: bool,
}

impl SlipBoundary {
    pub fn new(
        value: f64,
        dir: Direction,
        lower: bool,
        fixed: bool,
        no_slip: bool,
    ) -> SlipBoundary {
        SlipBoundary {
            line: BoundaryLine::<f64> { value, lower },
            direction: dir,
            fixed,
            no_slip,
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

#[derive(Debug, Clone)]
pub struct PeriodicBoundaryRect {
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

impl PeriodicBoundaryRect {
    pub fn new(x_min: f64, x_max: f64, y_min: f64, y_max: f64) -> Self {
        Self {
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }
}
