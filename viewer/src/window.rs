use std::sync::mpsc::{Receiver, Sender};

use mlsmpm::{Node, Particle};
use piston_window::*;

use crate::Snapshot;

pub struct ParticleWindow {
    width: u32,
    height: u32,

    space_size: f64,

    receiver: Receiver<Snapshot>,
    current: Snapshot,

    step_sender: Sender<usize>,

    previous_pressure_min: f64,
    previous_pressure_max: f64,
}

impl ParticleWindow {
    pub fn new(space_size: f64, receiver: Receiver<Snapshot>, step_sender: Sender<usize>) -> Self {
        ParticleWindow {
            width: 1280,
            height: 720,

            space_size,

            current: Snapshot::empty(),

            receiver,
            step_sender,

            previous_pressure_min: f64::MAX,
            previous_pressure_max: f64::MIN,
        }
    }

    pub fn run<F>(&mut self, data_seeker: &F)
    where
        F: Fn(&Snapshot),
    {
        Self::run_window(self.width, self.height, |event, window| {
            self.handle_step(event);
            self.draw_particles(event, window);
            self.snapshot(event, data_seeker);
            self.update_pressure_cache();
        });
    }

    fn run_window<F>(width: u32, height: u32, mut update: F)
    where
        F: FnMut(&Event, &mut PistonWindow),
    {
        let (width, height) = (width, height);
        let mut window: PistonWindow = WindowSettings::new("Hello Piston!", [width, height])
            .exit_on_esc(true)
            .build()
            .unwrap();
        while let Some(event) = window.next() {
            update(&event, &mut window);
        }
    }

    fn snapshot<F>(&self, event: &Event, data_seeker: &F)
    where
        F: Fn(&Snapshot),
    {
        if let Some(Button::Keyboard(key)) = event.press_args() {
            if key == Key::S {
                data_seeker(&self.current);
            }
        }
    }

    fn handle_step(&self, event: &Event) {
        if let Some(Button::Keyboard(key)) = event.press_args() {
            if key == Key::Space {
                self.step_sender.send(0).unwrap();
            } else if key == Key::Q {
                self.step_sender.send(1).unwrap();
            } else if key == Key::W {
                self.step_sender.send(10).unwrap();
            } else if key == Key::E {
                self.step_sender.send(100).unwrap();
            }
        }
    }

    fn draw_particles(&mut self, event: &Event, window: &mut PistonWindow) {
        let radius = 2.;
        let circle_rect = [-radius, -radius, radius, radius];

        let data_iter = self.receiver.try_iter();
        if let Some(snap) = data_iter.last() {
            self.current = snap;
        }

        window.draw_2d(event, |context, graphics, _device| {
            clear([1.0; 4], graphics);
            for p in self.current.particles.iter() {
                let v_p = convert_for_viewport(
                    (self.width, self.height),
                    self.space_size,
                    5.,
                    (p.x().x, p.x().y),
                );
                let transform = context.transform.trans(v_p.0, v_p.1);
                let color = convert_pressure_to_color(
                    p.v_norm(),
                    self.previous_pressure_min,
                    self.previous_pressure_max,
                );
                let color = [0., 0., 1., 1.];
                rectangle(color, circle_rect, transform, graphics);
            }
        });
    }

    fn update_pressure_cache(&mut self) {
        self.previous_pressure_min = f64::MAX;
        self.previous_pressure_max = f64::MIN;

        for p in self.current.particles.iter() {
            let pressure = p.v_norm();
            if pressure < self.previous_pressure_min {
                self.previous_pressure_min = pressure;
            }
            if pressure > self.previous_pressure_max {
                self.previous_pressure_max = pressure;
            }
        }
    }
}

fn convert_pressure_to_color(pressure: f64, min: f64, max: f64) -> [f32; 4] {
    let pressure = pressure as f32;
    let min = min as f32;
    let max = max as f32;
    let center = (max + min) / 2.;
    let r = f32::clamp((pressure - min) / (center - min), 0., 1.);
    let g = f32::clamp(1. - f32::abs((pressure - center) / (min - center)), 0., 1.);
    let b = f32::clamp((pressure - center) / (min - center), 0., 1.);

    [r, g, b, 1.]
}

fn convert_for_viewport(
    size: (u32, u32),
    space_size: f64,
    offset: f64,
    pos: (f64, f64),
) -> (f64, f64) {
    let (cx, cy) = (size.0 as f64 / 2., size.1 as f64 / 2.);
    let ratio = size.1 as f64 / space_size;

    (
        (pos.0 - offset) * ratio + cx,
        (pos.1 - offset) * ratio * -1. + cy,
    )
}
