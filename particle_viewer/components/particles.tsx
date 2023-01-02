import * as THREE from 'three'
import { createRoot } from 'react-dom/client'
import React, { useEffect, useMemo, useRef, useState } from 'react'
import { Canvas, useFrame, ThreeElements, useThree } from '@react-three/fiber'
import { OrthographicCamera } from '@react-three/drei'

export interface ParticleData {
	x: number[],
	v: number[],
}

interface ParticleProps {
	particles: ParticleData[]
}

function flattenPositions(plist: ParticleData[]): number[] {
	return plist.flatMap(p => [p.x[0] - 5, p.x[1] - 5, 0]);
}

export function Points() {
	const count = 100; // number point accross one axis ini akan generate point 10.00 dimana count hanya 100 karena multiply
	const sep = 3; //merupakan distance dari tiap point
	let positions = useMemo(() => {
		let positions = [];
		for (let xi = 0; xi < count; xi++) {
			for (let zi = 0; zi < count; zi++) {
				let x = sep * (xi - count / 2);
				let z = sep * (zi - count / 2);
				let y = 0;
				positions.push(x, y, z);
			}
		}
		return new Float32Array(positions); //merupakan array yang sesuai dengan buffer
	}, [count, sep]); //ini dibuat menjadi 1d array dikarenakan bufferAtribute tidak dapat menggunakan 2d array maka dari itu position array akan menjadi seperti [x1,y1,z1,x2,y2,z2,x....]
	return (
		<points>
			<bufferGeometry attach="geometry">
				<bufferAttribute
					attach="attributes-position" //attribute parameter yang akan dikontrol
					array={positions}
					count={positions.length / 3} //
					itemSize={3} //dikeranakan telah diketahui bahwa tiap arraytype axis akan berisi 3 value pada 1d array
				/>
			</bufferGeometry>
			<pointsMaterial
				attach="material"
				color={0x00aaff}
				size={0.5}
				sizeAttenuation //merupakan parameter yang menscale object berdasarkan perspective camera
				transparent={false}
				alphaTest={0.5} //merupakan thresshold saat rendering untuk mencega bila opacity dibawah value alphatest
				opacity={1.0}
			/>
		</points>
	);
}

export function Particles(props: ParticleProps) {
	const count = 100; // number point accross one axis ini akan generate point 10.00 dimana count hanya 100 karena multiply
	const sep = 3; //merupakan distance dari tiap point
	let dummy = useMemo(() => {
		let positions = [];
		for (let xi = 0; xi < count; xi++) {
			for (let zi = 0; zi < count; zi++) {
				let x = sep * (xi - count / 2);
				let z = sep * (zi - count / 2);
				let y = 0;
				positions.push(x, y, z);
			}
		}
		return new Float32Array(positions); //merupakan array yang sesuai dengan buffer
	}, [count, sep]);

	const positions = useMemo(() => {
		if (props.particles.length > 0) {
			const flatten = flattenPositions(props.particles);
			return new Float32Array(flatten);
		} else {
			return dummy;
		}
	}, [count, sep, props.particles]);

	// useThree(({ camera, size, gl }) => {
	// 	gl.setSize(size.width, size.height);
	// });

	const camera = useRef<THREE.OrthographicCamera>(null);

	useEffect(() => {
		if (camera == null || camera.current == null) {
			return;
		}

		camera.current.up.set(0, -1, 0);
		camera.current.position.set(0, 0, -10);
		camera.current.lookAt(new THREE.Vector3(0, 0, 0));
	}, [camera.current]);

	return (
		<OrthographicCamera ref={camera} makeDefault left={-1} right={1} top={-1} bottom={1}>
			<points>
				<bufferGeometry attach='geometry'>
					<bufferAttribute
						attach="attributes-position"
						array={positions}
						count={positions.length / 3}
						itemSize={3}
					/>
				</bufferGeometry>
				<pointsMaterial attach='material'
					color={0x00aaff}
					size={2} />
			</points>
		</OrthographicCamera>
	)
}

export function Box(props: ThreeElements['mesh']) {
	const mesh = useRef<THREE.Mesh>(null!)
	const [hovered, setHover] = useState(false)
	const [active, setActive] = useState(false)
	useFrame((state, delta) => (mesh.current.rotation.x += delta))
	return (
		<>
			<mesh
				{...props}
				ref={mesh}
				scale={active ? 1.5 : 1}
				onClick={(event) => setActive(!active)}
				onPointerOver={(event) => setHover(true)}
				onPointerOut={(event) => setHover(false)}>
				<boxGeometry args={[1, 1, 1]} />
				<meshStandardMaterial color={hovered ? 'hotpink' : 'orange'} />
			</mesh>
		</>
	)
}
