import { CATALOG, GROUPS } from './catalog-data.js';
import { convexHull, zonohedron, analyze } from './geom.js';
import { Viewer, solidObject, tilingObject, PALETTE } from './viewer.js';
import { generateTiles, checkCoverage, parallelohedronSpec } from './tiling.js';
import { setupPage, propsHTML, filterHTML, tilingHTML } from './ui.js';
import { t, L } from './i18n.js';

setupPage('catalog');

const $ = (s) => document.querySelector(s);
const side = $('#side'), stageTitle = $('#stage-title');
const viewer = new Viewer($('#stage'), { autoRotate: true });

const state = { mode: 'tiling', explode: 0.12, radius: 2.6, cut: 1, colorMode: 'mixed', params: {} };
let entry = null, current = null, queued = false, framed = false;

side.innerHTML = GROUPS.map((g) => `
  <h4>${L(g.title)}</h4><p class="gnote">${L(g.note)}</p>
  <nav>${CATALOG.filter((e) => e.group === g.id).map((e) =>
    `<a href="#${e.id}" data-id="${e.id}"><span class="dot ${e.status === 'no' ? 'no' : ''}"></span>${L(e.name)}</a>`).join('')}</nav>`).join('');

function specFor(e, poly, inf) {
  const s = e.tiling?.(state.params);
  if (!s) return null;
  if (s === 'parallelohedron') return parallelohedronSpec(poly, inf);
  return { kind: 'crystal', ...s };
}

function build() {
  queued = false;
  const e = entry;
  let poly, inf;
  try {
    poly = e.zonohedron ? zonohedron(e.zonohedron(state.params)) : convexHull(e.solid(state.params));
    inf = analyze(poly);
  } catch (err) {
    viewer.clear();
    current = null;
    $('#details').innerHTML = `<div class="verdict bad">${t(err.message)}</div>`;
    return;
  }
  // параметры вне семейства: оболочка уже не та (невыпуклое объединение) — разбиение не показываем
  const outside = e.expect && (!inf.equilateral || !e.expect.faces.includes(inf.nF));
  const spec = outside ? null : specFor(e, poly, inf);
  const mode = spec ? state.mode : 'solid';
  const R = state.radius * inf.circumradius;

  viewer.clear();
  let gen = null, cov = null;
  if (mode === 'tiling') {
    gen = generateTiles(inf, spec, R);
    current = tilingObject(poly, inf, gen);
    current.update(state);
    viewer.root.add(current.group);
    cov = checkCoverage(poly, inf, gen, R);
  } else {
    current = null;
    viewer.root.add(solidObject(poly, inf, { color: PALETTE[CATALOG.indexOf(e) % PALETTE.length] }));
  }
  if (!framed) { viewer.frame(mode === 'tiling' ? R + inf.circumradius : inf.circumradius * 1.25); framed = true; }

  document.querySelectorAll('.seg [data-mode]').forEach((b) => {
    b.classList.toggle('on', b.dataset.mode === mode);
    if (b.dataset.mode === 'tiling') b.disabled = !spec;
  });
  document.querySelectorAll('.tiling-only').forEach((x) => { x.style.display = mode === 'tiling' ? '' : 'none'; });

  const sub = e.status === 'no' ? t('st.notile') : t(mode === 'tiling' ? 'st.fragment' : 'st.solid');
  stageTitle.innerHTML = `<h2>${L(e.name)}</h2><div class="sub">${sub}</div>`;
  $('#details').innerHTML = `${outside ? `<div class="verdict bad" style="margin-bottom:12px">${t('warn.outside')}</div>` : ''}${propsHTML(inf)}${filterHTML(inf)}${mode === 'tiling' ? tilingHTML(gen, inf, cov) : ''}`;
}

function renderParams(e) {
  const params = (e.params || []).map((p) => `
    <div class="param"><label>${L(p.label)}<span class="mono" id="pv-${p.key}">${state.params[p.key]}${p.unit || ''}</span></label>
    <input type="range" min="${p.min}" max="${p.max}" step="${p.step}" value="${state.params[p.key]}" data-param="${p.key}"></div>`).join('');
  $('#about').innerHTML = `<p>${L(e.about)}</p>${params ? `<h3>${t('p.params')}</h3><div class="params">${params}</div>` : ''}`;
  $('#about').querySelectorAll('[data-param]').forEach((inp) => inp.addEventListener('input', () => {
    const p = e.params.find((q) => q.key === inp.dataset.param);
    state.params[p.key] = Number(inp.value);
    $(`#pv-${p.key}`).textContent = `${inp.value}${p.unit || ''}`;
    if (!queued) { queued = true; requestAnimationFrame(build); }
  }));
}

function select(id) {
  entry = CATALOG.find((e) => e.id === id) || CATALOG[0];
  state.params = Object.fromEntries((entry.params || []).map((p) => [p.key, p.value]));
  side.querySelectorAll('a').forEach((a) => a.classList.toggle('on', a.dataset.id === entry.id));
  side.querySelector('a.on')?.scrollIntoView({ block: 'nearest', inline: 'center' });
  framed = false;
  renderParams(entry);
  build();
}

document.querySelectorAll('.seg [data-mode]').forEach((b) => b.addEventListener('click', () => {
  state.mode = b.dataset.mode; framed = false; build();
}));
document.querySelectorAll('.seg [data-color]').forEach((b) => b.addEventListener('click', () => {
  state.colorMode = b.dataset.color;
  document.querySelectorAll('.seg [data-color]').forEach((x) => x.classList.toggle('on', x === b));
  current?.update(state);
}));
$('#explode').addEventListener('input', (ev) => { state.explode = Number(ev.target.value); current?.update(state); });
$('#cut').addEventListener('input', (ev) => { state.cut = Number(ev.target.value); current?.update(state); });
$('#radius').addEventListener('change', (ev) => { state.radius = Number(ev.target.value); framed = false; build(); });
$('#spin').addEventListener('change', (ev) => { viewer.controls.autoRotate = ev.target.checked; });

window.addEventListener('hashchange', () => select(location.hash.slice(1)));
select(location.hash.slice(1) || 'truncated-octahedron');
