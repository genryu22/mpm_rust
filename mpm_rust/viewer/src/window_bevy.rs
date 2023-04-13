use std::{
    f32::consts::{FRAC_PI_2, TAU},
    sync::mpsc,
};

use bevy::{
    prelude::*,
    render::{camera::ScalingMode, mesh},
    DefaultPlugins,
};
use bevy_points::{material::PointsShaderSettings, prelude::*};

use crate::Snapshot;

#[derive(Component, Default)]
struct Spiral;

#[derive(Component)]
struct Count(usize);

pub fn run(snapshot_receiver: mpsc::Receiver<Snapshot>) {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugin(PointsPlugin)
        .insert_resource(ClearColor(Color::rgb(0.01, 0.02, 0.08)))
        .insert_non_send_resource(snapshot_receiver)
        .add_startup_system(setup)
        //.add_system(test_animate)
        .add_system(update)
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<PointsMaterial>>,
) {
    commands.spawn(Count(0));

    commands.spawn((Spiral::default(), create_points_mesh(meshes, materials)));

    commands.spawn(Camera3dBundle {
        projection: OrthographicProjection {
            near: 0.1,
            far: 100.,
            scale: 1.0,
            scaling_mode: ScalingMode::WindowSize(180.),
            ..Default::default()
        }
        .into(),
        transform: Transform::from_translation(Vec3::new(0.0, 0.0, 5.0))
            .looking_at(Vec3::ZERO, Vec3::Y),
        ..Default::default()
    });
}

fn test_animate(
    mut meshes: ResMut<Assets<Mesh>>,
    mesh_query: Query<&Handle<Mesh>, With<Spiral>>,
    mut count: Query<&mut Count>,
) {
    if mesh_query.is_empty() {
        return;
    }
    let mesh_handle = mesh_query.single();
    if !meshes.contains(mesh_handle) {
        return;
    }
    let mesh = meshes.get_mut(mesh_handle).unwrap();
    let mut frame_count = count.get_single_mut().unwrap();
    let n = 320;
    let h = 3.0;
    update_mesh(
        mesh,
        PointsMesh::from_iter((0..n).map(|i| {
            let t01 = ((i as f32 + frame_count.0 as f32 / 10.) % n as f32) / ((n - 1) as f32);
            let r = t01 * TAU * 4.0;
            Vec3::new(r.cos(), (t01 - 0.5) * h, r.sin())
        })),
    );
    frame_count.0 += 1;
}

fn update_mesh(mesh: &mut Mesh, points_mesh: PointsMesh) {
    let vertices: Vec<[f32; 3]> = points_mesh
        .vertices
        .iter()
        .flat_map(|p| {
            let arr = p.to_array();
            [arr, arr, arr, arr]
        })
        .collect();
    let uv_set = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
    let uvs: Vec<[f32; 2]> = points_mesh.vertices.iter().flat_map(|_| uv_set).collect();
    let indices = bevy::render::mesh::Indices::U32(
        points_mesh
            .vertices
            .iter()
            .enumerate()
            .flat_map(|(i, _)| {
                let idx = (i * 4) as u32;
                [idx, idx + 1, idx + 3, idx + 2, idx + 3, idx + 1]
            })
            .collect(),
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, vertices);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
    if let Some(color) = points_mesh.colors {
        mesh.insert_attribute(
            Mesh::ATTRIBUTE_COLOR,
            color
                .iter()
                .flat_map(|c| {
                    let arr = c.as_rgba_f32();
                    [arr, arr, arr, arr]
                })
                .collect::<Vec<[f32; 4]>>(),
        );
    }
    mesh.set_indices(Some(indices));
}

fn update(world: &mut World) {
    let snapshot = {
        let receiver = world.get_non_send_resource::<mpsc::Receiver<Snapshot>>();
        if receiver.is_none() {
            return;
        }
        let received = receiver.unwrap().try_recv();

        received
    };
    if snapshot.is_err() {
        return;
    }
    let snapshot = snapshot.unwrap();

    let mesh_handle = world
        .query_filtered::<&Handle<Mesh>, With<Spiral>>()
        .get_single(world);
    if mesh_handle.is_err() {
        return;
    }
    let mesh_handle = mesh_handle.unwrap().clone();

    let meshes = world.get_resource_mut::<Assets<Mesh>>();
    if meshes.is_none() {
        return;
    }
    let mut meshes = meshes.unwrap();
    if !meshes.contains(&mesh_handle) {
        return;
    }
    let mesh = meshes.get_mut(&mesh_handle).unwrap();

    println!("{}", snapshot.steps);

    update_mesh(
        mesh,
        PointsMesh::from_iter(snapshot.particles.iter().map(|p| Vec3 {
            x: p.x().x as f32 - 3.,
            y: p.x().y as f32 - 3.,
            z: 0.,
        })),
    );
}

fn create_points_mesh(
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<PointsMaterial>>,
) -> MaterialMeshBundle<PointsMaterial> {
    let n = 320;
    let h = 3.0;

    MaterialMeshBundle {
        mesh: meshes.add(
            PointsMesh::from_iter((0..n).map(|i| {
                let t01 = (i as f32) / ((n - 1) as f32);
                let r = t01 * TAU * 4.0;
                Vec3::new(r.cos(), (t01 - 0.5) * h, r.sin())
            }))
            .into(),
        ),
        material: materials.add(PointsMaterial {
            settings: PointsShaderSettings {
                point_size: 20.,
                opacity: 1.,
                ..Default::default()
            },
            perspective: false,
            circle: true,
            ..Default::default()
        }),
        ..Default::default()
    }
}
