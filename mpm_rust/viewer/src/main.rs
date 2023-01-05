extern crate piston_window;

use std::{
    cell::RefCell,
    rc::Rc,
    sync::{mpsc, Arc, Mutex},
    thread,
};

use mlsmpm::Particle;
use piston_window::*;
use viewer::run_mlsmpm;

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
    let (tx, rx) = mpsc::channel();
    let handle = thread::spawn(move || {
        run_mlsmpm(tx);
    });

    let (width, height) = (640, 480);
    let mut window: PistonWindow = WindowSettings::new("Hello Piston!", [width, height])
        .exit_on_esc(true)
        .build()
        .unwrap();
    let particles: Rc<RefCell<Option<Vec<Particle>>>> = Rc::new(RefCell::new(Option::None));
    let rx_ref = Rc::new(rx);
    while let Some(event) = window.next() {
        let ps_ref = Rc::clone(&particles);
        let rx_ref = Rc::clone(&rx_ref);
        window.draw_2d(&event, move |context, graphics, _device| {
            clear([1.0; 4], graphics);
            if let Ok(ps) = rx_ref.try_recv() {
                *ps_ref.borrow_mut() = Some(ps);
            }
            if let Some(ps) = &*ps_ref.borrow() {
                for p in ps.iter() {
                    let v_p = convert_for_viewport((width, height), 1., 5., (p.x().x, p.x().y));
                    let radius = 10.;
                    ellipse(
                        [1.0, 0.0, 0.0, 1.0],
                        [v_p.0, v_p.1, radius, radius],
                        context.transform,
                        graphics,
                    );
                }
            }
        });
    }
}
