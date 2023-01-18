use std::error::Error;

use crate::Snapshot;

pub fn write_to_files(snapshot: &Snapshot) -> Result<(), Box<dyn Error>> {
    let now = chrono::Utc::now().format("%F-00-%M-%S").to_string();
    write_to_files_with_name(snapshot, &now)
}

pub fn write_to_files_with_name(snapshot: &Snapshot, name: &str) -> Result<(), Box<dyn Error>> {
    mlsmpm::file::write_particles(&snapshot.particles, snapshot.steps, "temp", name)?;
    mlsmpm::file::write_grid(&snapshot.grid, snapshot.steps, "temp", name)?;
    Ok(())
}
