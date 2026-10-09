import { zonohedron, analyze } from './geom.js';
import { Viewer, tilingObject } from './viewer.js';
import { generateTiles, parallelohedronSpec } from './tiling.js';
import { setupPage } from './ui.js';
import { L } from './i18n.js';
import { PROBLEMS, STATUS, PAGES } from './problems.js';

setupPage('home');

document.getElementById('problems').innerHTML = PROBLEMS.map(([title, desc, s, c], i) => `
  <div class="problem ${PAGES[i] ? 'active' : ''}"><div class="n">${String(i + 1).padStart(2, '0')}</div>
  <div><h3>${PAGES[i] ? `<a href="${PAGES[i]}">${L(title)}</a>` : L(title)}</h3><p>${L(desc)}</p></div><span class="chip ${c}">${L(STATUS[s])}</span></div>`).join('');

// Герой: дышащий фрагмент пены Кельвина из усечённых октаэдров
const viewer = new Viewer(document.getElementById('hero'), { autoRotate: true, interactive: false, fov: 30 });
viewer.controls.autoRotateSpeed = 0.35;
const gens = [[1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]].map((g) => g.map((x) => x / Math.SQRT2));
const poly = zonohedron(gens);
const inf = analyze(poly);
const R = 2.9 * inf.circumradius;
const obj = tilingObject(poly, inf, generateTiles(inf, parallelohedronSpec(poly, inf), R));
viewer.root.add(obj.group);
viewer.shiftX = (w) => (w > 900 ? 0.22 : 0);
viewer.resize();
viewer.frame(R * 1.3);
viewer.onFrame = (time) => obj.update({ explode: 0.16 + 0.12 * Math.sin(time / 1600), cut: 0.62 });
