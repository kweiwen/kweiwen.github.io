import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.118/build/three.module.js';
import {OrbitControls} from 'https://cdn.jsdelivr.net/npm/three@0.118/examples/jsm/controls/OrbitControls.js';
import {RectAreaLightHelper}  from 'https://cdn.jsdelivr.net/npm/three@0.118/examples/jsm/helpers/RectAreaLightHelper.js';
import {PositionalAudioHelper} from 'https://cdn.jsdelivr.net/npm/three@0.118/examples/jsm/helpers/PositionalAudioHelper.js';
import VS from './shaders/vertexShader.glsl.js';
import FS from './shaders/fragmentShader.glsl.js';

let DBG_CAM_POS = false;

let option = class {
  constructor() {
    this.message = 'dat.gui';
    this.rectLightWidth = 80;
    this.rectLightHeight = 80;
    this.rectLightIntensity = 0.85;
  }
};

class material {
  constructor() {

  }
}

class BasicWorld {
  constructor() {
    this.initAssets();
    this.initOptions();
    this.initRenderer();
    this.setCameraPositionAndControl();
    this.addLight();
    this.addTexture();
    this.addPlane();
    // let object1 = ;

    for(let col = 0; col < 7; col++)
    {
      for(let row = 0; row < 7; row++)
      {
        this.addListener(this.addSphere(4, 0, row*10 + 10, col*10 - 30));
      }
    }

    // let object2 = this.addSphere(4, 0, 30, 0);
    // this.addListener(object2);

    // let object3 = this.addSphere(4, 0, 40, 0);
    // this.addListener(object3);

    this.buffer = new Array(128).fill(0);

    window.addEventListener('resize', () => {
      this.OnWindowResize();
    }, false); 

    this.RAF();
  }

  initAssets() {
    this.previousRAF = null;
    this.totalTime = 0.0;
    this.index = 0;
    this.scene = new THREE.Scene();
    // const axesHelper = new THREE.AxesHelper(100);
    // this.scene.add(axesHelper);

    this._currentCameraPosition = new THREE.Vector3();
  }

  initOptions() {
    this.width = 10;
    this.height = 10;
  }

  initRenderer() {
    this.renderer = new THREE.WebGLRenderer({
      antialias: true,
    });

    this.renderer.shadowMap.enabled = true;
    this.renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(this.renderer.domElement);
  }

  setCameraPositionAndControl() {
    const fov = 45;
    const far = 1000.0;
    this.camera = new THREE.PerspectiveCamera(fov, window.innerWidth/window.innerHeight, 1, far);
    this.camera.position.set(150, 50, 0.0);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.target.set(0, 50, 0);
    this.controls.maxPolarAngle = THREE.Math.degToRad(100);
    this.controls.minPolarAngle = THREE.Math.degToRad(80);
    this.controls.maxDistance = 300;
    this.controls.minDistance = 50;
    this.controls.update();
    this.controls.addEventListener('change', this.onChangePosition.bind(this));
  }

  onChangePosition() {
    this.camera.getWorldPosition(this._currentCameraPosition);
    if(DBG_CAM_POS){console.log(this._currentCameraPosition);}
  }

  addLight() {
    const width = 80;
    const height = 80;
    const intensity = 0.85;
    this.rectLight1 = new THREE.RectAreaLight(0xffffff, intensity,  width, height);
    this.rectLight1.position.set(50, 40, -70);
    this.rectLight1.lookAt(0, 30, 0);
    this.scene.add(this.rectLight1);

    this.rectLight2 = new THREE.RectAreaLight(0xffffff, intensity,  width, height);
    this.rectLight2.position.set( 50, 40, 70);
    this.rectLight2.lookAt( 0, 30, 0 );
    this.scene.add(this.rectLight2);

    // this.rectLightHelper1 = new RectAreaLightHelper(this.rectLight1);
    // this.rectLight1.add(this.rectLightHelper1);
    // this.rectLightHelper2 = new RectAreaLightHelper(this.rectLight2);
    // this.rectLight2.add(this.rectLightHelper2);

  }

  addTexture() {
    const loader = new THREE.CubeTextureLoader();
    const texture = loader.load([
      // './resources/grid.jpg',
      // './resources/grid.jpg',
      // './resources/grid.jpg',
      // './resources/grid.jpg',
      // './resources/grid.jpg',
      // './resources/grid.jpg',
        // './resources/posx.jpg',
        // './resources/negx.jpg',
        // './resources/posy.jpg',
        // './resources/negy.jpg',
        // './resources/posz.jpg',
        // './resources/negz.jpg',
    ]);
    this.scene.background = texture;
  }

  addPlane() {
    const plane = new THREE.Mesh(
      new THREE.PlaneGeometry(512, 512, 0, 0),
      new THREE.MeshStandardMaterial({
        color: 0xFFFFFF,
      }));
    plane.castShadow = true;
    plane.receiveShadow = true;
    plane.rotation.x = -Math.PI / 2;
    plane.position.y = 0;
    this.scene.add(plane);
  }

  addSphere(radius, x, y, z) {
    let material = new THREE.ShaderMaterial({
      uniforms: {
        u_light_pos: {
          value: new THREE.Vector3(0, 0, 0),
        },
        u_noise_density: {
          value : 0.5,
        },
        u_noise_strength: {
          value : 0.5,
        },
        u_light_color: {
          value: new THREE.Color(0x888888),
        },
        u_light_intensity: {
          value : 0.65,
        },
        u_noise_coef: {
          value : 3.7,
        },
        u_noise_min: {
          value : 0.75,
        },
        u_noise_max: {
          value : 4,
        },
        u_noise_scale: {
          value : 0.8,
        },
        u_color: {
          value: new THREE.Color(0x666666),
        },
        u_time: {
          value: 0.0,
        },
        u_resolution: {
          value: new THREE.Vector2(window.innerWidth, window.innerHeight),
        },
      },
      vertexShader: VS,
      fragmentShader: FS,
    });

    let sphere = new THREE.Mesh(
      new THREE.SphereGeometry(radius, 64, 64),
      new THREE.ShaderMaterial(material)
      );
    sphere.material.uniforms.u_light_pos.value = new THREE.Vector3(x+30, y, z);
    // console.log(x+30, y, z);
    sphere.position.set(x, y, z);
    sphere.castShadow = true;
    sphere.receiveShadow = true;
    sphere.rotateY(THREE.Math.degToRad(90));
    this.scene.add(sphere);
    return sphere;
  }

  addListener(object) {
    this.listener = new THREE.AudioListener();
    this.camera.add(this.listener);
    this.listener.context.audioWorklet.addModule("digital-processor.js").then (() => this.addWorkletNode(object));
  }

  addWorkletNode(object) {
    const sound = new THREE.PositionalAudio(this.listener);
    let oscillator = new AudioWorkletNode(this.listener.context,  'drip-source-processor', { outputChannelCount: [1] });
    let level = new AudioWorkletNode(this.listener.context,  'level-meter-processor', { outputChannelCount: [1] });
    oscillator.connect(level);
    level.port.onmessage = (e) => {
      let energy = e.data.energy;
      let p2p = e.data.peak;
      // this.buffer.push(this.buffer.shift());
      // this.buffer[127] = energy * 64;
      // let peak = 0;
      // this.buffer.forEach((value, index, array) => {
      //   peak = peak + array[index];
      // })

      object.material.uniforms.u_noise_strength.value = Math.sqrt(Math.sqrt(Math.sqrt(p2p))) * 5;
      object.material.uniforms.u_noise_density.value = energy;    
      object.material.uniforms.u_time.value += 128 / 48000;
      object.material.uniforms.u_noise_scale.value = (Math.random() / 100) + 1.875;
    }
    sound.setNodeSource(level);
    sound.setRefDistance(10);
    sound.setDirectionalCone(180, 230, 0.1);
    // const helper = new PositionalAudioHelper(sound, 25);
    // sound.add(helper);
    object.add(sound);
  }  

  OnWindowResize() {
    this.camera.aspect = window.innerWidth / window.innerHeight;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(window.innerWidth, window.innerHeight);
  }

  RAF() {
    requestAnimationFrame((t) => {
      if (this.previousRAF === null) {
        this.previousRAF = t
      }

      this.RAF();
      this.renderer.render(this.scene, this.camera);
      this.Step(t - this.previousRAF);
      this.previousRAF = t;
    });
  }

  Step(timeElapsed) {
    const timeElapsedS = timeElapsed * 0.1;
    this.totalTime = this.totalTime + timeElapsedS;
    const tempTime = parseInt(this.totalTime);
  }
}

let APP = null;

window.addEventListener('DOMContentLoaded', () => {
  APP = new BasicWorld();
  // let gui = new dat.GUI({width:300});
  // let opts = new option();
  // let RectLightFolfer = gui.addFolder("RectLight");
  // let rectLightWidthSignal = RectLightFolfer.add(opts, 'rectLightWidth', 10, 100);
  // let rectLightHeightSignal = RectLightFolfer.add(opts, 'rectLightHeight', 10, 100);
  // let rectLightIntensitySignal = RectLightFolfer.add(opts, 'rectLightIntensity', 0, 1);

  // rectLightIntensitySignal.onChange(function(value) {
  //   APP.rectLight1.intensity = value;
  //   APP.rectLight2.intensity = value;
  //   APP.rectLightHelper1.update();
  //   APP.rectLightHelper2.update();
  // });

  // rectLightWidthSignal.onChange(function(value) {
  //   APP.rectLight1.width = value;
  //   APP.rectLight2.width = value;
  //   APP.rectLightHelper1.update();
  //   APP.rectLightHelper2.update();
  // });

  // rectLightHeightSignal.onChange(function(value) {
  //   APP.rectLight1.height = value;
  //   APP.rectLight2.height = value;
  //   APP.rectLightHelper1.update();
  //   APP.rectLightHelper2.update();
  // });
});