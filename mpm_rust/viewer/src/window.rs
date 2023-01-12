use std::sync::mpsc::{Receiver, Sender};

use mlsmpm::Particle;
use piston_window::*;

pub struct ParticleWindow {
    width: u32,
    height: u32,

    space_size: f64,

    receiver: Receiver<Vec<Particle>>,
    current_particles: Vec<Particle>,

    step_sender: Sender<usize>,
}

impl ParticleWindow {
    pub fn new(
        space_size: f64,
        receiver: Receiver<Vec<Particle>>,
        step_sender: Sender<usize>,
    ) -> Self {
        ParticleWindow {
            width: 1280,
            height: 720,

            space_size,

            current_particles: vec![],

            receiver,
            step_sender,
        }
    }

    pub fn run(&mut self) {
        Self::run_window(self.width, self.height, |event, window| {
            self.handle_step(event);
            self.draw_particles(event, window);
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
        let radius = 5.;
        let circle_rect = [-radius, -radius, radius, radius];
        let color = [1.0, 0.0, 0.0, 1.0];

        let ps_iter = self.receiver.try_iter();
        if let Some(particles) = ps_iter.last() {
            self.current_particles = particles;
        }

        window.draw_2d(event, |context, graphics, _device| {
            clear([1.0; 4], graphics);
            for p in self.current_particles.iter() {
                let v_p = convert_for_viewport(
                    (self.width, self.height),
                    self.space_size,
                    5.,
                    (p.x().x, p.x().y),
                );
                let transform = context.transform.trans(v_p.0, v_p.1);
                rectangle(color, circle_rect, transform, graphics);
            }
        });
    }
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
