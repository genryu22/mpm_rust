use mlsmpm::*;

fn main() {
    let settings = Settings {
        dt: 5e-3,
        gravity: -1e-2,
        dynamic_viscosity: 1e-2,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 200,
        rho_0: 1.,
        c: 0.,
        eos_power: 0.,
        boundary_mirror: false,
        vx_zero: true,
    };

    println!("{:?}", settings);

    let v_time_steps = (1. / settings.dynamic_viscosity / settings.dt).ceil() as u32;
    println!("粘性時間: L^2/mu = {} steps", v_time_steps);
    println!(
        "平均流速の理論解: U_mean = g*L*L/(nu*12) = {}",
        settings.gravity * 1. * 1. / ((settings.dynamic_viscosity / 1.) * 12.)
    );
    println!(
        "最大流速の理論解: U_max = 1.5*U_mean = {}",
        1.5 * settings.gravity * 1. * 1. / ((settings.dynamic_viscosity / 1.) * 12.)
    );

    let mut calc = Calculator::new(&settings, Space::new_for_poiseuille(&settings));
    calc.start(v_time_steps);

    println!("{}", calc.get_min_velocity());
}
