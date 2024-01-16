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
    pub effect_radius: usize,

    pub p2g_scheme: P2GSchemeType,
    pub g2p_scheme: G2PSchemeType,

    pub calc_convection_term: bool,

    pub pressure: Option<fn(&Particle, f64) -> f64>,
    pub pressure_grad: Option<fn(f64, f64, f64) -> Vector2<f64>>,

    pub reset_particle_position: bool,

    pub parallel: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            dt: 1e-4,
            gravity: 0.,
            dynamic_viscosity: 1e-2,
            alpha: 0.,
            affine: true,
            space_width: 10.,
            grid_width: 100,
            rho_0: 1.,
            c: 1e1,
            eos_power: 4.,
            boundary_mirror: false,
            vx_zero: false,
            weight_type: WeightType::QuadraticBSpline,
            effect_radius: 2,
            p2g_scheme: P2GSchemeType::MLSMPM,
            g2p_scheme: G2PSchemeType::MLSMPM,
            calc_convection_term: false,
            pressure: None,
            pressure_grad: None,
            reset_particle_position: false,
            parallel: true,
        }
    }
}

#[derive(Debug)]
pub enum WeightType {
    QuadraticBSpline,
    QuadraticBSpline2,
    CubicBSpline,
    CubicBSpline1_5,
    QuarticBSpline,
    QuinticBSpline,
    HexicBSpline,
    HepticBSpline,
    Linear,
    Spike,
}

#[derive(Debug, Clone, Copy)]
pub enum P2GSchemeType {
    MLSMPM,
    LSMPS,
    LsmpsLinear,
    Lsmps3rd,
    Lsmps4th,
    CompactLsmps,
    CompactLsmpsLinear,
    Compact0_1,
    Compact3_1,
    Compact0_2,
    Compact1_2,
    Compact2_2,
    Compact3_2,
    Compact3_3,
    Compact4_3,
    Compact4_4,
    CompactV0_1,
    CompactV1_1,
    CompactOnlyVelocity,
    CompactV3_1,
    CompactV0_2,
    CompactV1_2,
    CompactV2_2,
    CompactV3_2,
    CompactV0_3,
    CompactV1_3,
    CompactV2_3,
    CompactV3_3,
    CompactLaplacian2_2,
    CompactLaplacian3_2,
}

#[derive(Debug, Clone, Copy)]
pub enum G2PSchemeType {
    MLSMPM,
    LSMPS,
    Lsmps2ndMacro,
    Lsmps3rd,
    Lsmps4th,
    LsmpsLinear,
    LsmpsLinearMacro,
    CompactLsmps,
    LsmpsFlip1,
    LsmpsFlip2,
}

impl Settings {
    pub fn cell_width(&self) -> f64 {
        self.space_width / (self.grid_width as f64)
    }
}
