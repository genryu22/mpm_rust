const { createServer } = require('http')
const { parse } = require('url')
const next = require('next')

const PROTO_PATH = './protos/particle.proto';

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

const dev = process.env.NODE_ENV !== 'production'
const hostname = 'localhost'
const port = 3000
// when using middleware `hostname` and `port` must be provided below
const app = next({ dev, hostname, port })
const handle = app.getRequestHandler()

app.prepare().then(() => {
	const server = createServer(async (req, res) => {
		try {
			// Be sure to pass `true` as the second argument to `url.parse`.
			// This tells it to parse the query portion of the URL.
			const parsedUrl = parse(req.url, true)
			const { pathname, query } = parsedUrl

			await handle(req, res, parsedUrl)
		} catch (err) {
			console.error('Error occurred handling', req.url, err)
			res.statusCode = 500
			res.end('internal server error')
		}
	}).listen(port, (err) => {
		if (err) throw err
		console.log(`> Ready on http://${hostname}:${port}`)
	})

	const io = require('socket.io')(server);

	io.on('connection', (socket) => {
		console.log('A client connected.');
		socket.on('requestMPM', (payload) => {
			console.log(payload);

			const stream = client.RequestMPM({ name: 'test' });
			stream.on('data', (particles) => {
				socket.emit('particles-data', particles);
			})
		});
		socket.on('disconnect', () => {
			console.log('Conenction closed.');
		});
	});
})