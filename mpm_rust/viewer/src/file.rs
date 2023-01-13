use std::{
    error::Error,
    fs::{self, File},
    path::Path,
};

use crate::Snapshot;

pub fn write_to_files(snapshot: &Snapshot) -> Result<(), Box<dyn Error>> {
    let folder = Path::new("temp");
    if !folder.exists() {
        fs::create_dir(folder)?;
    }

    let now = chrono::Utc::now().format("%F-00-%M-%S").to_string();

    let mut grid =
        csv::Writer::from_path(folder.join(format!("{}-{}-grid.csv", now, snapshot.steps)))?;
    grid.write_record(&[
        "i", "vx", "vy", "v_starx", "v_stary", "forcex", "forcey", "mass",
    ])?;
    for (i, n) in snapshot.grid.iter().enumerate() {
        let formatted = n.formatted_list();
        grid.write_record(&[
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
    grid.flush()?;

    let mut particle =
        csv::Writer::from_path(folder.join(format!("{}-{}-particle.csv", now, snapshot.steps)))?;
    particle.write_record(&["x", "y", "vx", "vy", "mass"])?;
    for p in snapshot.particles.iter() {
        let formatted = p.formatted_list();
        particle.write_record(&[
            &formatted[0],
            &formatted[1],
            &formatted[2],
            &formatted[3],
            &formatted[4],
        ])?;
    }
    grid.flush()?;

    Ok(())
}
