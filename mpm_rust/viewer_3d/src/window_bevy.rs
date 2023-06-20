use std::{
    f32::consts::{FRAC_PI_2, TAU},
    sync::mpsc,
};

use bevy::{
    input::mouse::MouseMotion,
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

pub fn run(snapshot_receiver: mpsc::Receiver<Snapshot>, store_data: fn(&Snapshot)) {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                resolution: bevy::window::WindowResolution::new(1920., 1080.),
                ..default()
            }),
            ..default()
        }))
        .add_plugin(PointsPlugin)
        .insert_resource(ClearColor(Color::rgb(0.01, 0.02, 0.08)))
        .insert_non_send_resource(snapshot_receiver)
        .insert_non_send_resource(store_data)
        .add_startup_system(setup)
        .add_system(update)
        .add_system(move_camera)
        .add_system(rotate_camera)
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
        projection: PerspectiveProjection {
            near: 0.1,
            far: 100.,
            ..Default::default()
        }
        .into(),
        transform: Transform::from_translation(Vec3::new(0.0, 0.0, 5.0))
            .looking_at(Vec3::ZERO, Vec3::Y),
        ..Default::default()
    });

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
}

fn move_camera(
    mut camera_transform: Query<&mut Transform, With<Camera>>,
    keys: Res<Input<KeyCode>>,
) {
    if camera_transform.is_empty() {
        return;
    }
    let mut camera_transform = camera_transform.single_mut();

    {
        let mut delta = Vec3::ZERO;
        if keys.pressed(KeyCode::W) {
            delta += camera_transform.forward();
        }
        if keys.pressed(KeyCode::A) {
            delta += camera_transform.left();
        }
        if keys.pressed(KeyCode::S) {
            delta += camera_transform.back();
        }
        if keys.pressed(KeyCode::D) {
            delta += camera_transform.right();
        }
        if keys.pressed(KeyCode::LShift) {
            delta += camera_transform.down();
        }
        if keys.pressed(KeyCode::Space) {
            delta += camera_transform.up();
        }
        camera_transform.translation += 0.1 * delta;
    }
}

fn rotate_camera(
    mut camera_transform: Query<&mut Transform, With<Camera>>,
    buttons: Res<Input<MouseButton>>,
    mut motion_evr: EventReader<MouseMotion>,
    primary_window: Query<&Window, With<bevy::window::PrimaryWindow>>,
) {
    if camera_transform.is_empty() || primary_window.is_empty() {
        return;
    }
    let mut camera_transform = camera_transform.single_mut();
    let window = primary_window.single();

    if !buttons.pressed(MouseButton::Left) {
        return;
    }

    let (mut yaw, mut pitch, _) = camera_transform.rotation.to_euler(EulerRot::YXZ);

    for ev in motion_evr.iter() {
        let sens = 0.002;
        let window_scale = window.height().min(window.width());
        pitch -= (sens * ev.delta.y * window_scale).to_radians();
        yaw -= (sens * ev.delta.x * window_scale).to_radians();
    }

    pitch = pitch.clamp(-1.54, 1.54);

    // Order is important to prevent unintended roll
    camera_transform.rotation =
        Quat::from_axis_angle(Vec3::Y, yaw) * Quat::from_axis_angle(Vec3::X, pitch);
}

fn update(world: &mut World) {
    let snapshot = {
        let receiver = world.get_non_send_resource::<mpsc::Receiver<Snapshot>>();
        if receiver.is_none() {
            return;
        }

        if false {
            receiver.unwrap().try_recv()
        } else {
            receiver
                .unwrap()
                .try_iter()
                .last()
                .map_or(Err(mpsc::TryRecvError::Empty), |s| Ok(s))
        }
    };

    if snapshot.is_err() {
        return;
    }
    let snapshot = snapshot.unwrap();

    let store_data = world.get_non_send_resource::<fn(&Snapshot)>();
    if let Some(store_data) = store_data {
        store_data(&snapshot);
    }

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

    let mut points_mesh = PointsMesh::from_iter(snapshot.particles.iter().map(|p| Vec3 {
        x: p.x().x as f32 - 5.,
        y: p.x().y as f32 - 5.,
        z: p.x().z as f32 - 5.,
    }));

    fn convert_particle_to_scaler(particle: &mlsmpm_3d::Particle) -> f64 {
        particle.pressure()
    }

    let (scaler_min, scaler_max) = {
        let mut min = f64::MAX;
        let mut max = f64::MIN;

        for p in snapshot.particles.iter() {
            min = f64::min(min, convert_particle_to_scaler(p));
            max = f64::max(max, convert_particle_to_scaler(p));
        }

        (min, max)
    };

    points_mesh.colors = Some(
        snapshot
            .particles
            .iter()
            .map(|p| {
                let color =
                    convert_scalar_to_color(convert_particle_to_scaler(p), scaler_min, scaler_max);

                Color::Rgba {
                    red: color[0],
                    green: color[1],
                    blue: color[2],
                    alpha: color[3],
                }
            })
            .collect(),
    );

    update_mesh(mesh, points_mesh);

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

    fn convert_scalar_to_color(scalar: f64, min: f64, max: f64) -> [f32; 4] {
        let pressure = scalar as f32;
        let min = min as f32;
        let max = max as f32;
        let center = (max + min) / 2.;
        let r = f32::clamp((pressure - min) / (center - min), 0., 1.);
        let g = f32::clamp(1. - f32::abs((pressure - center) / (min - center)), 0., 1.);
        let b = f32::clamp((pressure - center) / (min - center), 0., 1.);

        [r, g, b, 1.]
    }
}