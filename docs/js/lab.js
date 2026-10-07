import { convexHull, zonohedron, analyze, mcmullen, mergeGenerators, v } from './geom.js';
import { Viewer, solidObject, tilingObject, PALETTE } from './viewer.js';
import { generateTiles, checkCoverage, parallelohedronSpec } from './tiling.js';
import { setupPage, propsHTML, filterHTML, tilingHTML, mcmullenHTML } from './ui.js';
import { t, L } from './i18n.js';

setupPage('lab');

const $ = (s) => document.querySelector(s);
const viewer = new Viewer($('#stage'));
const phi = (1 + Math.sqrt(5)) / 2;
const h = Math.sqrt(3) / 2;
const fmt = (rows) => rows.map((r) => r.map((x) => +x.toFixed(6)).join(' ')).join('\n');

const ZONO = [
  [{ en: 'Cube', ru: 'Куб', es: 'Cubo' }, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]],
  [{ en: 'Rhombohedron', ru: 'Ромбоэдр', es: 'Romboedro' }, [[1, 0, 0], [0.5, h, 0], [0.5, h / 3, Math.sqrt(2 / 3)]]],
  [{ en: 'Hex prism', ru: 'Шестиуг. призма', es: 'Prisma hex.' }, [[1, 0, 0], [0.5, h, 0], [-0.5, h, 0], [0, 0, 1]]],
  [{ en: 'Rhombic dodecahedron', ru: 'Ромбододекаэдр', es: 'Dodecaedro rómbico' }, [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]],
  [{ en: 'Bilinski', ru: 'Билинский', es: 'Bilinski' }, [[0, 1, phi], [0, -1, phi], [1, phi, 0], [-1, phi, 0]]],
  [{ en: 'Elongated dodecahedron', ru: 'Удл. додекаэдр', es: 'Dodecaedro alargado' }, [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1], [0, 0, Math.sqrt(3)]]],
  [{ en: 'Truncated octahedron', ru: 'Усеч. октаэдр', es: 'Octaedro truncado' }, [[1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]]],
  [{ en: 'Rhombic icosahedron', ru: 'Ромбоикосаэдр', es: 'Icosaedro rómbico' }, [[0, 1, phi], [0, -1, phi], [1, phi, 0], [-1, phi, 0], [phi, 0, 1]]],
  [{ en: 'Triacontahedron', ru: 'Триаконтаэдр', es: 'Triacontaedro' }, [[0, 1, phi], [0, -1, phi], [1, phi, 0], [-1, phi, 0], [phi, 0, 1], [phi, 0, -1]]],
];

const pm = [-1, 1];
const HULL = [
  [{ en: 'Tetrahedron', ru: 'Тетраэдр', es: 'Tetraedro' }, [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]],
  [{ en: 'Octahedron', ru: 'Октаэдр', es: 'Octaedro' }, [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]],
  [{ en: 'Icosahedron', ru: 'Икосаэдр', es: 'Icosaedro' }, pm.flatMap((a) => pm.flatMap((b) => [[0, a, b * phi], [a, b * phi, 0], [a * phi, 0, b]]))],
  [{ en: 'Dodecahedron', ru: 'Додекаэдр', es: 'Dodecaedro' }, [
    ...pm.flatMap((a) => pm.flatMap((b) => pm.map((c) => [a, b, c]))),
    ...pm.flatMap((a) => pm.flatMap((b) => [[0, a / phi, b * phi], [a / phi, b * phi, 0], [a * phi, 0, b / phi]])),
  ]],
  [{ en: 'Pyramid J1', ru: 'Пирамида J1', es: 'Pirámide J1' }, [[0.5, 0.5, 0], [-0.5, 0.5, 0], [-0.5, -0.5, 0], [0.5, -0.5, 0], [0, 0, Math.SQRT1_2]]],
  [{ en: 'Gyrobifastigium', ru: 'Гиробифастигиум', es: 'Girobifastigio' }, [[0.5, 0.5, 0], [-0.5, 0.5, 0], [-0.5, -0.5, 0], [0.5, -0.5, 0], [0.5, 0, h], [-0.5, 0, h], [0, 0.5, -h], [0, -0.5, -h]]],
  [{ en: 'Cuboctahedron', ru: 'Кубооктаэдр', es: 'Cuboctaedro' }, pm.flatMap((a) => pm.flatMap((b) => [[a, b, 0], [a, 0, b], [0, a, b]]))],
];

const state = { tool: 'zono', mode: 'solid', explode: 0.12, cut: 1 };
let current = null;

function parse(text) {
  const rows = text.split('\n').map((l) => l.trim()).filter((l) => l && !l.startsWith('#'))
    .map((l) => l.split(/[\s,;]+/).map(Number));
  if (rows.some((r) => r.length !== 3 || r.some((x) => !Number.isFinite(x)))) throw new Error('err.parse');
  return rows;
}

function presets(el, list, input, run) {
  el.innerHTML = list.map(([name], i) => `<button data-i="${i}">${L(name)}</button>`).join('');
  el.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => { input.value = fmt(list[b.dataset.i][1]); run(); }));
}

function show(poly, title, extraHTML, spec) {
  const inf = analyze(poly);
  const mode = spec ? state.mode : 'solid';
  viewer.clear();
  let gen = null, cov = null;
  if (mode === 'tiling') {
    const R = 2.6 * inf.circumradius;
    gen = generateTiles(inf, spec, R);
    current = tilingObject(poly, inf, gen);
    current.update(state);
    viewer.root.add(current.group);
    cov = checkCoverage(poly, inf, gen, R);
    viewer.frame(R + inf.circumradius);
  } else {
    current = null;
    viewer.root.add(solidObject(poly, inf, { color: PALETTE[state.tool === 'zono' ? 0 : 1] }));
    viewer.frame(inf.circumradius * 1.25);
  }
  document.querySelectorAll('.toolbar [data-mode]').forEach((b) => {
    b.classList.toggle('on', b.dataset.mode === mode);
    if (b.dataset.mode === 'tiling') b.disabled = !spec;
  });
  document.querySelectorAll('.tiling-only').forEach((x) => { x.style.display = mode === 'tiling' ? '' : 'none'; });
  $('#stage-title').innerHTML = `<h2>${title}</h2><div class="sub">${t('st.faces', { n: inf.nF })} · ${t(inf.equilateral ? 'st.equi' : 'st.nonequi')}</div>`;
  $('#info').innerHTML = `${extraHTML || ''}${propsHTML(inf)}${filterHTML(inf)}${mode === 'tiling' ? tilingHTML(gen, inf, cov) : ''}`;
}

const fail = (e) => { $('#err').textContent = t(e.message); };

function runZono() {
  try {
    $('#err').textContent = '';
    let gens = parse($('#zono-input').value);
    if ($('#zono-norm').checked) gens = gens.map((g) => (v.len(g) ? v.unit(g) : g));
    const merged = mergeGenerators(gens);
    const poly = zonohedron(merged);
    const mc = mcmullen(merged);
    const spec = mc.ok ? parallelohedronSpec(poly, analyze(poly)) : null;
    show(poly, t('lab.zono'), mcmullenHTML(mc, merged.length), spec);
  } catch (e) { fail(e); }
}

function runHull() {
  try {
    $('#err').textContent = '';
    show(convexHull(parse($('#hull-input').value)), t('lab.title.hull'), '', null);
  } catch (e) { fail(e); }
}

presets($('#zono-presets'), ZONO, $('#zono-input'), runZono);
presets($('#hull-presets'), HULL, $('#hull-input'), runHull);
$('#zono-go').addEventListener('click', runZono);
$('#hull-go').addEventListener('click', runHull);
$('#zono-rand').addEventListener('click', () => {
  const r = () => { let p; do { p = [0, 0, 0].map(() => Math.random() * 2 - 1); } while (v.len(p) > 1 || v.len(p) < 0.2); return v.unit(p); };
  $('#zono-input').value = fmt([r(), r(), r(), r()]);
  runZono();
});

document.querySelectorAll('.tabs [data-tool]').forEach((b) => b.addEventListener('click', () => {
  state.tool = b.dataset.tool;
  state.mode = 'solid';
  document.querySelectorAll('.tabs [data-tool]').forEach((x) => x.classList.toggle('on', x === b));
  document.querySelectorAll('[data-panel]').forEach((p) => { p.hidden = p.dataset.panel !== state.tool; });
  (state.tool === 'zono' ? runZono : runHull)();
}));
document.querySelectorAll('.toolbar [data-mode]').forEach((b) => b.addEventListener('click', () => {
  state.mode = b.dataset.mode;
  (state.tool === 'zono' ? runZono : runHull)();
}));
$('#explode').addEventListener('input', (e) => { state.explode = Number(e.target.value); current?.update(state); });
$('#cut').addEventListener('input', (e) => { state.cut = Number(e.target.value); current?.update(state); });
$('#spin').addEventListener('change', (e) => { viewer.controls.autoRotate = e.target.checked; });

$('#zono-input').value = fmt(ZONO[4][1]);
$('#hull-input').value = fmt(HULL[5][1]);
runZono();
