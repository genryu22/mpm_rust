use std::{
    error::Error,
    fs::{self, File},
    path::Path,
};

use crate::{Node, Particle};

pub fn write_particles(
    particles: &Vec<Particle>,
    steps: usize,
    folder_name: &str,
    name: &str,
) -> Result<(), Box<dyn Error>> {
    let folder = Path::new(folder_name);
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let mut particle_writer =
        csv::Writer::from_path(folder.join(format!("{}-{}-particle.csv", name, steps)))?;
    particle_writer.write_record(&["x", "y", "vx", "vy", "mass", "pressure"])?;
    for p in particles.iter() {
        let formatted = p.formatted_list();
        particle_writer.write_record(&[
            &formatted[0],
            &formatted[1],
            &formatted[2],
            &formatted[3],
            &formatted[4],
            &formatted[5],
        ])?;
    }
    particle_writer.flush()?;

    Ok(())
}

pub fn write_grid(
    grid: &Vec<Node>,
    steps: usize,
    folder_name: &str,
    name: &str,
) -> Result<(), Box<dyn Error>> {
    let folder = Path::new(folder_name);
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let mut grid_writer =
        csv::Writer::from_path(folder.join(format!("{}-{}-grid.csv", name, steps)))?;
    grid_writer.write_record(&[
        "i", "vx", "vy", "v_starx", "v_stary", "forcex", "forcey", "mass",
    ])?;
    for (i, n) in grid.iter().enumerate() {
        let formatted = n.formatted_list();
        grid_writer.write_record(&[
            &i.to_string(),
            &formatted[0],
            &formatted[1],
            &formatted[2],
            &formatted[3],
            &formatted[4],
            &formatted[5],
            &formatted[6],
        ])?;
    }
    grid_writer.flush()?;

    Ok(())
}
