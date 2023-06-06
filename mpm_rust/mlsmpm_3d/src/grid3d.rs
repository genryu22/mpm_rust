use std::{
    cell::RefCell,
    vec::{self, IntoIter},
};

use crate::*;

#[derive(Debug)]
pub struct Grid {
    pub nodes: Vec<Vec<Vec<Node>>>,
    pub slip_walls: Vec<SlipWall>,

    pub size: usize,
}

impl Grid {
    pub fn get_node(&self, index: Vector3i) -> Option<&Node> {
        if index.x < 0 || index.y < 0 || index.z < 0 {
            return None;
        }

        self.nodes.get(index.z as usize).and_then(|yx| {
            yx.get(index.y as usize)
                .and_then(|x| x.get(index.x as usize))
        })
    }
    pub fn get_node_mut(&mut self, index: Vector3i) -> Option<&mut Node> {
        if index.x < 0 || index.y < 0 || index.z < 0 {
            return None;
        }

        self.nodes.get_mut(index.z as usize).and_then(|yx| {
            yx.get_mut(index.y as usize)
                .and_then(|x| x.get_mut(index.x as usize))
        })
    }

    pub fn all_nodes(&self) -> Vec<&Node> {
        self.all_nodes_indices()
            .iter()
            .map(|(x, y, z)| &self.nodes[*z][*y][*x])
            .collect::<Vec<_>>()
    }

    pub fn all_nodes_indices(&self) -> Vec<(usize, usize, usize)> {
        (0..self.size)
            .flat_map(|x| (0..self.size).flat_map(move |y| (0..self.size).map(move |z| (x, y, z))))
            .collect::<Vec<_>>()
    }

    pub fn all_nodes_iter(&mut self) -> GridMutIterator {
        GridMutIterator::new(self)
    }
}
pub struct GridMutIterator<'a> {
    grid: &'a mut Grid,
    indices_iter: IntoIter<(usize, usize, usize)>,
}

impl<'a> Iterator for GridMutIterator<'a> {
    type Item = &'a mut Node;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((x, y, z)) = self.indices_iter.next() {
            Some(unsafe {
                self.grid
                    .nodes
                    .as_mut_ptr()
                    .add(z)
                    .as_mut()
                    .unwrap()
                    .as_mut_ptr()
                    .add(y)
                    .as_mut()
                    .unwrap()
                    .as_mut_ptr()
                    .add(x)
                    .as_mut()
                    .unwrap()
            })
        } else {
            None
        }
    }
}

impl<'a> GridMutIterator<'a> {
    pub fn new(grid: &'a mut Grid) -> GridMutIterator<'a> {
        Self {
            indices_iter: grid.all_nodes_indices().into_iter(),
            grid,
        }
    }
}

#[derive(Debug)]
pub struct SlipWall {
    pub center: Vector3f,
    pub a: Vector3f,
    pub b: Vector3f,
    pub n: Vector3f,
}

impl SlipWall {
    pub fn new(center: Vector3f, a: Vector3f, b: Vector3f) -> Self {
        Self {
            center,
            a,
            b,
            n: calc_n(&center, &a, &b),
        }
    }

    pub fn calc_distance(&self, target: &Vector3f) -> f64 {
        self.n.dot(&(target - self.center))
    }

    pub fn calc_projection(&self, target: &Vector3f) -> Vector3f {
        target - self.calc_distance(target) * self.n
    }
}

fn calc_n(center: &Vector3f, a: &Vector3f, b: &Vector3f) -> Vector3f {
    (a - center).cross(&(b - center)).normalize()
}
