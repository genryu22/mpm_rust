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

export function ParticlesRaw(props: ParticleProps) {
	const rootRef = useRef<HTMLDivElement>(null);
	const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
	const cameraRef = useRef<THREE.OrthographicCamera | null>(null);

	useEffect(() => {
		const renderer = new THREE.WebGLRenderer();
		rendererRef.current = renderer;
		if (rootRef.current != null) {
			rootRef.current.appendChild(renderer.domElement);
		}

		renderer.setSize(window.innerWidth, window.innerHeight);

		const SIZE = 1;
		let camera_width;
		let camera_height;
		if (window.innerWidth <= window.innerHeight) {
			camera_width = SIZE;
			camera_height = SIZE / window.innerWidth * window.innerHeight;
		} else {
			camera_width = SIZE / window.innerHeight * window.innerWidth;
			camera_height = SIZE;
		}
		const camera = new THREE.OrthographicCamera(camera_width / -2, camera_width / 2, camera_height / -2, camera_height / 2);
		camera.up.set(0, -1, 0);
		camera.position.set(0, 0, -10);
		camera.lookAt(new THREE.Vector3(0, 0, 0));
		cameraRef.current = camera;

		return () => {
			if (rootRef.current != null) {
				rootRef.current.removeChild(renderer.domElement);
				rendererRef.current = null;
				cameraRef.current = null;
			}
		}
	}, [rootRef.current]);

	const sceneRef = useRef<THREE.Scene | null>(null);

	useEffect(() => {
		const scene = new THREE.Scene();
		sceneRef.current = scene;
		const geometry = new THREE.BufferGeometry();
		geometry.setAttribute('position', new THREE.Float32BufferAttribute(flattenPositions(props.particles), 3));

		const material = new THREE.PointsMaterial({
			size: 2,
			color: 0xffffff,
		});

		const pointsMesh = new THREE.Points(geometry, material);
		scene.add(pointsMesh);
	}, [props.particles]);

	useEffect(() => {
		let canceled = false;
		function animate() {
			if (canceled) {
				return;
			}
			rendererRef.current?.render(sceneRef.current!, cameraRef.current!);
			requestAnimationFrame(animate);
		}
		animate();

		return () => {
			canceled = true;
		}
	});

	return <div ref={rootRef}></div>
}