use std::{
    f32::consts::{FRAC_PI_2, TAU},
    sync::mpsc,
};

use bevy::{prelude::*, DefaultPlugins};
use bevy_points::{material::PointsShaderSettings, prelude::*};

use crate::Snapshot;

pub fn run(snapshot_receiver: mpsc::Receiver<Snapshot>) {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugin(PointsPlugin)
        .insert_resource(ClearColor(Color::rgb(0.01, 0.02, 0.08)))
        .insert_non_send_resource(snapshot_receiver)
        .add_startup_system(setup)
        .add_system(update)
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<PointsMaterial>>,
) {
    let ORIGIN: Vec3 = Vec3::new(0.0, 0.0, -5.0);
    let mut pt = PointsMesh::from(Mesh::from(shape::UVSphere {
        radius: 1.0,
        sectors: 36,
        stacks: 18,
    }));
    let n = pt.vertices.len() as f32;
    pt.colors = Some(
        pt.vertices
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let t = (i as f32) / (n - 1.) * 360.;
                Color::hsl(t, 1., 0.5)
            })
            .collect(),
    );

    commands.spawn(MaterialMeshBundle {
        mesh: meshes.add(pt.into()),
        material: materials.add(PointsMaterial {
            settings: PointsShaderSettings {
                point_size: 0.1,
                opacity: 0.5,
                ..Default::default()
            },
            perspective: true,
            alpha_mode: AlphaMode::Blend,
            ..Default::default()
        }),
        transform: Transform::from_translation(Vec3::NEG_X * 1.25)
            .with_rotation(Quat::from_axis_angle(Vec3::ONE.normalize(), FRAC_PI_2)),
        ..Default::default()
    });

    let n = 320;
    let h = 3.0;
    commands.spawn(MaterialMeshBundle {
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
        transform: Transform::from_translation(Vec3::X * 1.25),
        ..Default::default()
    });

    commands.spawn(Camera3dBundle {
        projection: PerspectiveProjection {
            fov: 45.0,
            aspect_ratio: 1.,
            near: 0.1,
            far: 100.,
        }
        .into(),
        transform: Transform::from_translation(ORIGIN).looking_at(Vec3::ZERO, Vec3::Y),
        ..Default::default()
    });
}

fn update(world: &mut World) {
    let receiver = world.get_non_send_resource::<mpsc::Receiver<Snapshot>>();
    if receiver.is_none() {
        return;
    }
    let receiver = receiver.unwrap();

    if let Ok(snapshot) = receiver.try_recv() {
        println!("{}", snapshot.steps);
    }
}
