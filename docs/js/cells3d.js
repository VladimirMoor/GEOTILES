// 3D-поверхности настоящих клеток (пузыри, зёрна, клетки тканей, скутоиды) из data/cells/*.json.
// Разметка: <div class="cells3d" data-src="data/cells/scutoid.json"></div> — одна сцена с несколькими клетками;
// <div class="gallery cellgal" data-src="data/cells/gallery.json"></div> — по карточке на каждую клетку.
// JSON: {"meshes":[{"v":[x,y,z,…],"f":[i,j,k,…],"color":"#…","name":{en,ru,es},"note":{…},"opacity":1}]};
// координаты уже отцентрированы и нормированы на радиус ≈ 1. Вьюер создаётся, когда блок становится видимым.
import * as THREE from 'three';
import { Viewer } from './viewer.js';
import { L } from './i18n.js';

function meshObject(m) {
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(m.v, 3));
  g.setIndex(m.f);
  if (m.c) g.setAttribute('color', new THREE.Float32BufferAttribute(m.c, 3));
  g.computeVertexNormals();
  const op = m.opacity ?? 1;
  const mat = new THREE.MeshStandardMaterial({
    color: m.c ? '#ffffff' : (m.color || '#4f8fd9'), vertexColors: !!m.c, roughness: 0.55, metalness: 0.05, flatShading: !!m.flat,
    transparent: op < 1, opacity: op, side: THREE.DoubleSide, depthWrite: op >= 1,
  });
  const obj = new THREE.Group();
  obj.add(new THREE.Mesh(g, mat));
  if (m.flat) obj.add(new THREE.LineSegments(new THREE.EdgesGeometry(g, 20), new THREE.LineBasicMaterial({ color: '#111' })));
  if (m.label) obj.add(textSprite(L(m.label), m.labelPos || [0, 0, 0]));
  return obj;
}

// подпись, всегда видимая поверх модели
function textSprite(text, pos) {
  const c = document.createElement('canvas');
  const ctx = c.getContext('2d');
  const font = 'bold 44px Manrope, sans-serif';
  ctx.font = font;
  c.width = Math.ceil(ctx.measureText(text).width) + 24; c.height = 64;
  ctx.font = font; ctx.textBaseline = 'middle';
  ctx.lineWidth = 8; ctx.strokeStyle = 'rgba(0,0,0,0.75)'; ctx.strokeText(text, 12, 33);
  ctx.fillStyle = '#ffffff'; ctx.fillText(text, 12, 33);
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, depthTest: false, transparent: true }));
  s.renderOrder = 10;
  const h = 0.16;
  s.scale.set(h * c.width / c.height, h, 1);
  s.position.set(...pos);
  return s;
}

function radiusOf(meshes) {
  let r = 0;
  const solid = meshes.filter((m) => (m.opacity ?? 1) >= 1);  // кадр по самим клеткам, без вспомогательных плоскостей
  for (const m of solid.length ? solid : meshes) for (let i = 0; i < m.v.length; i += 3) r = Math.max(r, Math.hypot(m.v[i], m.v[i + 1], m.v[i + 2]));
  return r || 1;
}

function scene(stage, meshes, { interactive = true, speed = 0.8 } = {}) {
  const v = new Viewer(stage, { autoRotate: true, interactive });
  v.controls.autoRotateSpeed = speed;
  for (const m of meshes) v.root.add(meshObject(m));
  v.frame(radiusOf(meshes) * 1.15);
  return v;
}

function whenVisible(el, fn) {
  const io = new IntersectionObserver((es) => {
    for (const e of es) if (e.isIntersecting) { io.disconnect(); fn(); }
  }, { rootMargin: '200px' });
  io.observe(el);
}

export function setupCells3d() {
  for (const el of document.querySelectorAll('.cells3d[data-src]')) {
    whenVisible(el, async () => {
      const d = await (await fetch(el.dataset.src)).json();
      const stage = document.createElement('div');
      stage.className = 'cells3d-stage';
      el.prepend(stage);
      scene(stage, d.meshes);
    });
  }
  for (const el of document.querySelectorAll('.cellgal[data-src]')) {
    whenVisible(el, async () => {
      const d = await (await fetch(el.dataset.src)).json();
      for (const m of d.meshes) {
        const card = document.createElement('div');
        card.className = 'mini';
        card.innerHTML = `<div class="mini-stage"></div><span class="mini-cap"><b>${L(m.name)}</b>${m.note ? `<span>${L(m.note)}</span>` : ''}</span>`;
        el.appendChild(card);
        scene(card.querySelector('.mini-stage'), [m], { speed: 1.6 });
      }
    });
  }
}
