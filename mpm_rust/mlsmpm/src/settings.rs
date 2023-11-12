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

    pub pressure: Option<fn(&Particle, f64) -> f64>,

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
            pressure: None,
            reset_particle_position: false,
            parallel: true,
        }
    }
}

#[derive(Debug)]
pub enum WeightType {
    QuadraticBSpline,
    CubicBSpline,
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
    Compact_0_1,
    Compact_3_1,
    Compact_0_2,
    Compact_1_2,
    Compact_2_2,
    Compact_3_2,
    Compact_3_3,
    Compact_4_3,
    Compact_4_4,
    Compact_v_0_1,
    Compact_v_1_1,
    CompactOnlyVelocity,
    Compact_v_3_1,
    Compact_v_0_2,
    Compact_v_1_2,
    Compact_v_2_2,
    Compact_v_3_2,
    Compact_v_0_3,
    Compact_v_1_3,
    Compact_v_2_3,
    Compact_v_3_3,
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
}

impl Settings {
    pub fn cell_width(&self) -> f64 {
        self.space_width / (self.grid_width as f64)
    }
}
