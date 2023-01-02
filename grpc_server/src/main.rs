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
        Err(Status::unavailable(""))
    }
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let addr = "[::1]:50051".parse()?;
    // let greeter = MyGreeter::default();

    // Server::builder()
    //     .add_service(GreeterServer::new(greeter))
    //     .serve(addr)
    //     .await?;

    Ok(())
}
