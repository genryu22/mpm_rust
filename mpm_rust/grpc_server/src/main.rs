use mlsmpm::{Calculator, Particle, Settings, Space};
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tonic::{transport::Server, Request, Response, Status};

use particle::mpm_streamer_server::{MpmStreamer, MpmStreamerServer};
use particle::{MpmRequest, Particles};

pub mod particle {
    tonic::include_proto!("particle"); // The string specified here must match the proto package name
}

#[derive(Debug, Default)]
pub struct MyParticleStreamer {}

#[tonic::async_trait]
impl MpmStreamer for MyParticleStreamer {
    type RequestMPMStream = ReceiverStream<Result<Particles, Status>>;
    async fn request_mpm(
        &self,
        request: Request<MpmRequest>,
    ) -> Result<Response<Self::RequestMPMStream>, Status> {
        let (tx, rx) = mpsc::channel(128);
        tokio::spawn(async move {
            let settings = Settings {
                dt: 0.05,
                gravity: 1e-2,
                dynamic_viscosity: 1e-2,
                alpha: 0.,
                affine: true,
                space_width: 10.,
                grid_width: 200,
            };

            println!("{:?}", settings);

            let mut calc = Calculator::new(&settings, Space::new_for_poiseuille(&settings));

            for _i in 0..100 {
                calc.update();
                tx.send(Ok(create_particles_packet(calc.get_particles())))
                    .await
                    .unwrap();
            }
        });

        Ok(Response::new(ReceiverStream::new(rx)))
    }
}

fn create_particles_packet(particles: &Vec<Particle>) -> Particles {
    let mut x = Vec::<f64>::with_capacity(particles.len() * 2);
    let mut y = Vec::<f64>::with_capacity(particles.len() * 2);
    let mut vx = Vec::<f64>::with_capacity(particles.len() * 2);
    let mut vy = Vec::<f64>::with_capacity(particles.len() * 2);

    for p in particles.iter() {
        let pos = p.x();
        let vel = p.v();
        x.push(pos.x);
        y.push(pos.y);
        vx.push(vel.x);
        vy.push(vel.y);
    }

    Particles { x, y, vx, vy }
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let addr = "[::1]:50051".parse()?;
    let greeter = MyParticleStreamer::default();

    Server::builder()
        .add_service(MpmStreamerServer::new(greeter))
        .serve(addr)
        .await?;

    Ok(())
}
