use std::{net::TcpListener, thread::spawn};

use tungstenite::{
    accept, accept_hdr,
    handshake::{server::Request, server::Response},
};

pub struct WSServer {}

impl WSServer {
    pub fn start() {
        let server = TcpListener::bind("0.0.0.0:9001").unwrap();
        for stream in server.incoming() {
            spawn(move || {
                let websocket = accept(stream.unwrap());
                if let Ok(mut websocket) = websocket {
                    loop {
                        let msg = websocket.read_message().unwrap();

                        // We do not want to send back ping/pong messages.
                        if msg.is_binary() || msg.is_text() {
                            websocket.write_message(msg).unwrap();
                        }
                    }
                }
            });
        }
    }
}
