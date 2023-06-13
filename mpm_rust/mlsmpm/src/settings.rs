use crate::*;

#[derive(Debug)]
pub struct Settings {
    pub dt: f64,
    pub gravity: f64,
    pub dynamic_viscosity: f64,
    pub alpha: f64,
    pub affine: bool,
    pub space_width: f64,
    pub grid_width: U,

    pub rho_0: f64,

    pub c: f64,
    pub eos_power: f64,

    pub boundary_mirror: bool,
    pub vx_zero: bool,

    pub weight_type: WeightType,
}

#[derive(Debug)]
pub enum WeightType {
    QuadraticBSpline,
    QubicBSpline,
    Linear,
}

impl Settings {
    pub fn cell_width(&self) -> f64 {
        self.space_width / (self.grid_width as f64)
    }
}
