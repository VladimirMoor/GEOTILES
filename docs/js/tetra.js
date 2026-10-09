// Страница задачи 2: 40 спорадических рациональных тетраэдров и рёбра без звезды.
import * as THREE from 'three';
import { TETRA } from './tetra-data.js';
import { convexHull, analyze } from './geom.js';
import { Viewer, solidObject, PALETTE } from './viewer.js';
import { setupPage, esc } from './ui.js';
import { L } from './i18n.js';

setupPage('research');

const T = {
  pick: { en: 'Pick a tetrahedron', ru: 'Выберите тетраэдр', es: 'Elige un tetraedro' },
  angles: { en: 'Dihedral angles (×π)', ru: 'Двугранные углы (×π)', es: 'Ángulos diedros (×π)' },
  edges: { en: 'Edges (unit volume)', ru: 'Рёбра (единичный объём)', es: 'Aristas (volumen unitario)' },
  bad: { en: 'edge without a face-to-face star', ru: 'ребро без звезды «лицом к лицу»', es: 'arista sin estrella cara a cara' },
  status: { en: 'Status', ru: 'Статус', es: 'Estado' },
  tiles: { en: 'tiles (classical)', ru: 'замощает (классический)', es: 'tesela (clásico)' },
  nof2f: { en: 'no face-to-face tiling; non-face-to-face: open', ru: 'нет разбиения «лицом к лицу»; не «лицом к лицу» — открыто', es: 'sin teselación cara a cara; no cara a cara: abierto' },
  isnew: { en: 'new (not covered by Chentouf–Sun, Prop. 2.7)', ru: 'новое (не покрыто предл. 2.7 Chentouf–Sun)', es: 'nuevo (no cubierto por la Prop. 2.7 de Chentouf–Sun)' },
  old: { en: 'already follows from Chentouf–Sun, Prop. 2.7', ru: 'следовало из предл. 2.7 Chentouf–Sun', es: 'ya se sigue de la Prop. 2.7 de Chentouf–Sun' },
  clusters: { en: 'Convex clusters of ≤ 7 copies', ru: 'Выпуклых кластеров из ≤ 7 копий', es: 'Clústeres convexos de ≤ 7 copias' },
  legendBad: { en: 'red: edges without a star', ru: 'красные — рёбра без звезды', es: 'rojo: aristas sin estrella' },
};

const $ = (s) => document.querySelector(s);
const viewer = new Viewer($('#tstage'), { autoRotate: true });
const EDGES = ['12', '34', '13', '24', '14', '23'];

$('#tpick').innerHTML = `<div class="small muted" style="margin-bottom:8px">${L(T.pick)}</div>` + TETRA.map((r) =>
  `<button class="tbtn ${r.status.startsWith('Sommerville') ? 'ok' : r.f2f_exclusion_new ? 'new' : ''}" data-no="${r.no}">#${r.no}</button>`).join('');

function cylinder(a, b, rad, color) {
  const A = new THREE.Vector3(...a), B = new THREE.Vector3(...b);
  const len = A.distanceTo(B);
  const g = new THREE.CylinderGeometry(rad, rad, len, 10);
  const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ color }));
  m.position.copy(A.clone().add(B).multiplyScalar(0.5));
  m.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), B.clone().sub(A).normalize());
  return m;
}

function show(no) {
  const r = TETRA.find((x) => x.no === no);
  document.querySelectorAll('.tbtn').forEach((b) => b.classList.toggle('on', +b.dataset.no === no));
  viewer.clear();
  const poly = convexHull(r.V);
  const inf = analyze(poly);
  const obj = solidObject(poly, inf, { color: r.status.startsWith('Sommerville') ? PALETTE[2] : PALETTE[0], opacity: 0.55 });
  for (const e of r.edges_without_star) {
    const i = +e[0] - 1, j = +e[1] - 1;
    obj.add(cylinder(r.V[i], r.V[j], inf.minLen * 0.035, '#d8605d'));
  }
  viewer.root.add(obj);
  viewer.frame(inf.circumradius * 1.25);
  const known = r.status.startsWith('Sommerville');
  $('#tinfo').innerHTML = `
    <h3>#${r.no}</h3>
    <dl class="kv">
      <dt>${L(T.status)}</dt><dd><span class="chip ${known ? 'ok' : 'info'}">${known ? esc(r.status.replace(' (tiles)', '')) + ' — ' + L(T.tiles) : L(T.nof2f)}</span></dd>
      ${known ? '' : `<dt></dt><dd class="small">${L(r.f2f_exclusion_new ? T.isnew : T.old)}</dd>`}
      <dt>${L(T.clusters)}</dt><dd>${r.convex_clusters_le7 ?? '—'}</dd>
    </dl>
    <h3>${L(T.angles)}</h3>
    <div class="small mono">(α12, α34, α13, α24, α14, α23) = (${r.angles.join(', ')})</div>
    <h3>${L(T.edges)}</h3>
    <dl class="kv small">${EDGES.map((e) => `<dt>${e}</dt><dd>${r.len[e]}${r.edges_without_star.includes(e) ? ` <span class="chip bad">${L(T.bad)}</span>` : ''}</dd>`).join('')}</dl>
    <p class="small muted">${L(T.legendBad)}</p>`;
}

$('#tpick').addEventListener('click', (ev) => {
  const b = ev.target.closest('.tbtn');
  if (b) show(+b.dataset.no);
});
show(25);
