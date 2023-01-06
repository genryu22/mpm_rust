extern crate piston_window;

use std::{
    cell::RefCell,
    rc::Rc,
    sync::{mpsc, Arc, Mutex},
    thread,
};

use mlsmpm::Particle;
use piston_window::*;
use viewer::{run_dambreak, run_mlsmpm};

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

fn main() {
    let particles: Arc<Mutex<Vec<Particle>>> = Arc::new(Mutex::new(vec![]));

    let (tx, rx) = mpsc::channel();

    let rx_op = Some(rx);
    let rx_op = None;

    {
        let particles = Arc::clone(&particles);
        let handle = thread::spawn(move || {
            //run_dambreak(particles, Some(rx));
            run_mlsmpm(particles, rx_op);
        });
    }

    let (width, height) = (640, 480);
    let mut window: PistonWindow = WindowSettings::new("Hello Piston!", [width, height])
        .exit_on_esc(true)
        .build()
        .unwrap();
    let radius = 5.;
    let circle_rect = [-radius, -radius, radius, radius];
    let color = [1.0, 0.0, 0.0, 1.0];
    while let Some(event) = window.next() {
        let particles = Arc::clone(&particles);
        let tx = mpsc::Sender::clone(&tx);

        if let Some(Button::Keyboard(key)) = event.press_args() {
            if key == Key::Space {
                tx.send(true).unwrap();
            }
        }

        window.draw_2d(&event, move |context, graphics, _device| {
            clear([1.0; 4], graphics);
            if let Ok(ref mut particles) = particles.try_lock() {
                for p in particles.iter() {
                    let v_p = convert_for_viewport((width, height), 1., 5., (p.x().x, p.x().y));
                    let transform = context.transform.trans(v_p.0, v_p.1);
                    rectangle(color, circle_rect, transform, graphics);
                }
            }
        });
    }
}
