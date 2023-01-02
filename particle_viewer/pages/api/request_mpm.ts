// Next.js API route support: https://nextjs.org/docs/api-routes/introduction
import type { NextApiRequest, NextApiResponse } from 'next'
import path from 'path';

const PROTO_PATH = path.join(process.cwd(), './protos/particle.proto');

var parseArgs = require('minimist');
var grpc = require('@grpc/grpc-js');
var protoLoader = require('@grpc/proto-loader');
var packageDefinition = protoLoader.loadSync(
	PROTO_PATH,
	{
		keepCase: true,
		longs: String,
		enums: String,
		defaults: true,
		oneofs: true
	});
var particle_proto = grpc.loadPackageDefinition(packageDefinition).particle;

const target = 'localhost:50051';
const client = new particle_proto.MPMStreamer(target, grpc.credentials.createInsecure());

type Data = {
	name: string
}

export default function handler(
	req: NextApiRequest,
	res: NextApiResponse<Data>
) {
	const stream = client.RequestMPM({ name: 'test' });
	stream.on('data', (particles: any) => {
		console.log(particles);
	})

	res.status(200).json({ name: 'John Doe' })
}