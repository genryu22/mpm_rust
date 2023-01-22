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
}

impl Settings {
    pub fn new(
        dt: f64,
        gravity: f64,
        dynamic_viscosity: f64,
        alpha: f64,
        affine: bool,
        space_width: f64,
        grid_width: U,
        rho_0: f64,
        c: f64,
        eos_power: f64,
        boundary_mirror: bool,
        vx_zero: bool,
    ) -> Self {
        Self {
            dt,
            gravity,
            dynamic_viscosity,
            alpha,
            affine,
            space_width,
            grid_width,
            rho_0,
            c,
            eos_power,
            boundary_mirror,
            vx_zero,
        }
    }

    pub fn cell_width(&self) -> f64 {
        self.space_width / (self.grid_width as f64)
    }
}
