use mlsmpm::{Calculator, Particle, Settings, Space};
use tokio::sync::mpsc;
use tokio_stream::wrappers::ReceiverStream;
use tonic::{transport::Server, Request, Response, Status};

use particle::mpm_streamer_server::{MpmStreamer, MpmStreamerServer};
use particle::{MpmRequest, Particles};

pub mod particle {
    tonic::include_proto!("particle"); // The string specified here must match the proto package name

    pub(crate) const FILE_DESCRIPTOR_SET: &[u8] =
        tonic::include_file_descriptor_set!("particle_descriptor");
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
        let mpm_request = request.into_inner();
        let (tx, rx) = mpsc::channel(1);
        tokio::spawn(async move {
            let settings = Settings {
                dt: mpm_request.dt,
                gravity: mpm_request.gravity,
                dynamic_viscosity: mpm_request.dynamic_viscosity,
                alpha: mpm_request.alpha,
                affine: mpm_request.affine,
                space_width: mpm_request.space_width,
                grid_width: mpm_request.grid_width as usize,
            };

            println!("{:?}", settings);

            let mut calc = Calculator::new(&settings, Space::new_for_poiseuille(&settings));

            for _i in 0..mpm_request.step_count {
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

    let service = tonic_reflection::server::Builder::configure()
        .register_encoded_file_descriptor_set(particle::FILE_DESCRIPTOR_SET)
        .build()
        .unwrap();

    Server::builder()
        .add_service(service)
        .add_service(MpmStreamerServer::new(greeter))
        .serve(addr)
        .await?;

    Ok(())
}
