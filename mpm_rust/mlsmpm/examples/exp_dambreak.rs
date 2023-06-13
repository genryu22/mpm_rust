use std::{error::Error, fs, path::Path, thread, vec};

use mlsmpm::*;

fn main() -> Result<(), Box<dyn Error>> {
    let settings = Settings {
        dt: 1e-4,
        gravity: -10.,
        dynamic_viscosity: 1e-3,
        alpha: 0.,
        affine: true,
        space_width: 10.,
        grid_width: 500,
        rho_0: 1000.,
        c: 1e1,
        eos_power: 4.,
        boundary_mirror: true,
        vx_zero: false,
        weight_type: WeightType::QuadraticBSpline,
    };

    println!("{:?}", settings);

    let folder = Path::new("exp_dambreak");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let mut particle_writer = csv::Writer::from_path(folder.join("dambreak.csv"))?;
    particle_writer.write_record(&["nt", "z"])?;

    let mut calc = Calculator::new(&settings, Space::new_for_dambreak_experiment(&settings));

    particle_writer.write_record(&["0", &(calc.get_max_x() - 3.).to_string()])?;

    for i in 0..7900 {
        calc.update();

        let step = (i + 1) as f64;
        particle_writer.write_record(&[
            &(step * settings.dt * f64::sqrt(-2. * settings.gravity)).to_string(),
            &(calc.get_max_x() - 3.).to_string(),
        ])?;
    }

    particle_writer.flush()?;

    Ok(())
}
