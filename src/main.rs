use mpm_rust::{Calculator, Settings, Space};

fn main() {
    let settings = Settings {
        dt: 0.05,
        gravity: 1e-2,
        dynamic_viscosity: 1e-2,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 200,
    };

    println!("{:?}", settings);

    let mut calc = Calculator::new(&settings, Space::new_for_poiseuille(&settings));
    calc.start(100);

    //println!("{:?}", calc);
}
