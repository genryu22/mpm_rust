use nalgebra::{Matrix3, Vector3};

type U = usize;

type Vector3f = Vector3<f64>;
type Vector3i = Vector3<i64>;
type Vector3u = Vector3<U>;
type Matrix3f = Matrix3<f64>;

mod bspline_iter;
mod grid3d;
mod node3d;
mod particle3d;
mod settings;
mod space3d;

pub use bspline_iter::*;
pub use grid3d::*;
pub use node3d::*;
pub use particle3d::*;
pub use settings::*;
