// three.js-вьюер: одна сцена на контейнер, математическая ось z смотрит вверх.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

export const PALETTE = ['#4f8fd9', '#e3a33b', '#4fb38a', '#d8605d', '#8d78d8', '#3fb2c4', '#d983b0', '#9aaa45'];

const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

export class Viewer {
  constructor(el, { autoRotate = true, interactive = true, fov = 35 } = {}) {
    this.el = el;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    el.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(fov, 1, 0.01, 1000);
    this.camera.position.set(4.2, 2.6, 5.2);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.autoRotate = autoRotate;
    this.controls.autoRotateSpeed = 0.7;
    this.controls.enabled = interactive;

    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x3a4552, 1.25));
    const key = new THREE.DirectionalLight(0xffffff, 1.7);
    key.position.set(5, 8, 6);
    const rim = new THREE.DirectionalLight(0xffffff, 0.55);
    rim.position.set(-6, -3, -5);
    this.scene.add(key, rim);

    this.root = new THREE.Group();
    this.root.rotation.x = -Math.PI / 2;
    this.scene.add(this.root);

    this.onFrame = null;
    new ResizeObserver(() => this.resize()).observe(el);
    this.resize();
    const loop = (t) => {
      requestAnimationFrame(loop);
      this.controls.update();
      this.onFrame?.(t);
      this.renderer.render(this.scene, this.camera);
    };
    requestAnimationFrame(loop);
  }

  resize() {
    const w = this.el.clientWidth || 1, h = this.el.clientHeight || 1;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    // shiftX: сдвиг картинки вправо в долях ширины (на главной слева текст)
    const s = typeof this.shiftX === 'function' ? this.shiftX(w) : 0;
    if (s) this.camera.setViewOffset(w, h, -w * s, 0, w, h); else this.camera.clearViewOffset();
    this.camera.updateProjectionMatrix();
  }

  clear() {
    for (const ch of [...this.root.children]) {
      this.root.remove(ch);
      ch.traverse((o) => { o.geometry?.dispose?.(); });
    }
  }

  frame(radius) {
    const fov = (this.camera.fov * Math.PI) / 180;
    const dist = (radius / Math.sin(fov / 2)) * 1.08 / Math.min(1, this.camera.aspect);
    const dir = this.camera.position.clone().sub(this.controls.target).normalize();
    this.controls.target.set(0, 0, 0);
    this.camera.position.copy(dir.multiplyScalar(dist));
    this.camera.near = dist / 100;
    this.camera.far = dist * 20;
    this.camera.updateProjectionMatrix();
  }
}

export function polyGeometry(poly) {
  const pos = [], nor = [];
  poly.F.forEach((f, fi) => {
    const n = poly.N[fi];
    for (let i = 1; i < f.length - 1; i++) for (const k of [f[0], f[i], f[i + 1]]) { pos.push(...poly.V[k]); nor.push(...n); }
  });
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(nor, 3));
  return g;
}

export function edgeGeometry(poly, info) {
  const pos = [];
  for (const e of info.edges) pos.push(...poly.V[e.a], ...poly.V[e.b]);
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  return g;
}

const faceMaterial = (color, opacity = 1) => new THREE.MeshStandardMaterial({
  color, flatShading: true, roughness: 0.5, metalness: 0.05,
  transparent: opacity < 1, opacity, side: opacity < 1 ? THREE.DoubleSide : THREE.FrontSide,
  polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1,
});

/** Одно тело: грани, рёбра и вершины. */
export function solidObject(poly, info, { color = PALETTE[0], opacity = 0.92 } = {}) {
  const group = new THREE.Group();
  group.add(new THREE.Mesh(polyGeometry(poly), faceMaterial(color, opacity)));
  group.add(new THREE.LineSegments(edgeGeometry(poly, info), new THREE.LineBasicMaterial({ color: cssVar('--edge') || '#111' })));
  const r = Math.max(0.012, info.minLen * 0.035);
  const sg = new THREE.SphereGeometry(r, 12, 8);
  const sm = new THREE.MeshStandardMaterial({ color: cssVar('--edge') || '#111' });
  for (const p of poly.V) { const s = new THREE.Mesh(sg, sm); s.position.set(...p); group.add(s); }
  group.position.set(...info.centroid.map((x) => -x));
  return group;
}

/** Фрагмент разбиения. Возвращает группу и функцию update({explode, cut, colorMode}). */
export function tilingObject(poly, info, gen) {
  const group = new THREE.Group();
  const geo = polyGeometry(poly), egeo = edgeGeometry(poly, info);
  const mats = PALETTE.map((c) => faceMaterial(c));
  const emat = new THREE.LineBasicMaterial({ color: cssVar('--edge') || '#111', transparent: true, opacity: 0.55 });
  const items = gen.tiles.map((t) => {
    const mesh = new THREE.Mesh(geo, mats[t.cls]);
    const lines = new THREE.LineSegments(egeo, emat);
    mesh.matrixAutoUpdate = lines.matrixAutoUpdate = false;
    group.add(mesh, lines);
    return { t, mesh, lines, off: t.centroid.clone().sub(gen.c0) };
  });
  group.position.copy(gen.c0).multiplyScalar(-1);
  const zs = items.map((it) => it.off.z);
  const zmin = Math.min(...zs), zmax = Math.max(...zs);
  const shift = new THREE.Matrix4();
  const update = ({ explode = 0, cut = 1, colorMode = 'mixed' } = {}) => {
    const zcut = zmin + (zmax - zmin) * cut + 1e-6;
    for (const it of items) {
      shift.makeTranslation(it.off.x * explode, it.off.y * explode, it.off.z * explode).multiply(it.t.matrix);
      it.mesh.matrix.copy(shift);
      it.lines.matrix.copy(shift);
      it.mesh.visible = it.lines.visible = it.off.z <= zcut;
      it.mesh.material = mats[colorMode === 'orbit' ? it.t.orbit % mats.length : it.t.cls];
    }
  };
  update();
  return { group, update, count: items.length };
}
