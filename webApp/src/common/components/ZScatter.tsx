import { range } from 'd3';
import React from 'react';
import { compose } from 'recompose';
import {
  Camera,
  CanvasTexture,
  Geometry,
  Line,
  LineBasicMaterial,
  PerspectiveCamera,
  Scene,
  Sprite,
  SpriteMaterial,
  Vector3,
  WebGLRenderer,
} from 'three';
import { OrbitControls } from 'three-orbitcontrols-ts';

const PI2 = Math.PI * 2;

class ZScatter extends React.Component<any, any> {
  private mount: HTMLDivElement;
  private scene: Scene;
  private camera: Camera;
  private renderer: WebGLRenderer;
  private frameRequestId: number;
  private frame: number;
  private width: number;
  private height: number;

  constructor(props) {
    super(props);
    this.frame = 0;
    this.width = props.width;
    this.height = props.height;

    this.renderer = new WebGLRenderer({ alpha: true, antialias: true });
    this.renderer.setClearColor(0x000000, 0);
    this.renderer.setSize(this.width, this.height);
  }

  public render() {
    const { className } = this.props;
    return (
      <div
        style={{ width: this.width, height: this.height }}
        className={className}
        ref={(mount) => {
          this.mount = mount;
        }}
      />
    );
  }

  public componentDidMount() {
    this.plot();

    this.mount.appendChild(this.renderer.domElement);
    this.start();
  }

  public componentDidUpdate() {
    /**
     * Note that if React wants to reuse the same ZScatter instance, we need to reset the scene, but we can
     * still reuse the renderer.
     *
     * This also means that if any prop changes (xs, ys, zs, width, or height) the entire scene is rebuilt.
     * Later we could probably be smarter about it.
     *
     */
    this.plot();
  }

  public componentWillUnmount() {
    this.stop();
    this.mount.removeChild(this.renderer.domElement);
  }

  private plot = () => {
    this.width = this.mount.clientWidth;
    this.height = this.mount.clientHeight;

    const canvasDot: HTMLCanvasElement = document.createElement('canvas');
    canvasDot.height = 128;
    canvasDot.width = 128;
    const contextDot = canvasDot.getContext('2d');
    contextDot.lineWidth = 3;
    contextDot.strokeStyle = 'gray';
    contextDot.fillStyle = 'white';
    contextDot.beginPath();
    contextDot.arc(64, 64, 60, 0, PI2, true);
    contextDot.fill();
    contextDot.stroke();

    this.scene = new Scene();
    range(5000).forEach(() => {
      const materialDot = new SpriteMaterial({
        map: new CanvasTexture(canvasDot),
        color: Math.random() * 0x808008 + 0x808080,
      });
      const spriteDot = new Sprite(materialDot);
      spriteDot.position.x = Math.random() * 800 - 400;
      spriteDot.position.y = Math.random() * 800 - 400;
      spriteDot.position.z = Math.random() * 800 - 400;
      spriteDot.scale.x = spriteDot.scale.y = 10;
      this.scene.add(spriteDot);
    });

    const materialLine = new LineBasicMaterial({ color: 0x777777 });
    const geomFrameLines = new Geometry();
    geomFrameLines.vertices.push(new Vector3(-400, -400, -400));
    geomFrameLines.vertices.push(new Vector3(-400, -400, 400));
    geomFrameLines.vertices.push(new Vector3(400, -400, 400));
    geomFrameLines.vertices.push(new Vector3(400, -400, -400));
    geomFrameLines.vertices.push(new Vector3(-400, -400, -400));
    const line = new Line(geomFrameLines, materialLine);
    this.scene.add(line);

    this.camera = new PerspectiveCamera(75, this.width / this.height, 0.1, 5000);
    this.camera.position.set(1000, 0, 0);

    const controls = new OrbitControls(this.camera, this.renderer.domElement);
    controls.enableZoom = true;
    controls.target = new Vector3(0, 0, 0);
  };

  private start = () => {
    if (!this.frameRequestId) {
      this.frameRequestId = window.requestAnimationFrame(this.animate);
    }
  };

  private stop = () => {
    window.cancelAnimationFrame(this.frameRequestId);
  };

  private animate = () => {
    // this.camera.position.x = Math.cos(this.frame / 500) * 1000;
    // this.camera.position.z = Math.sin(this.frame / 500) * 1000;
    // this.frame += 1;
    // this.camera.lookAt(new Vector3(0, 0, 0));

    this.renderScene();
    this.frameRequestId = window.requestAnimationFrame(this.animate);
  };

  private renderScene = () => {
    this.renderer.render(this.scene, this.camera);
  };
}

const enhance: any = compose();

export default enhance(ZScatter);
