// Маленькие вращающиеся 3D-модели из каталога внутри текстовых страниц.
// Разметка: <div class="gallery"><a class="mini" data-id="engel38"></a>…</div> — подпись и ссылка на каталог
// добавляются сами. Вьюер создаётся, только когда карточка попадает в зону видимости.
import { CATALOG } from './catalog-data.js';
import { convexHull, zonohedron, analyze } from './geom.js';
import { Viewer, solidObject, PALETTE } from './viewer.js';
import { L } from './i18n.js';

const FACES = { en: 'faces', ru: 'граней', es: 'caras' };

function polyOf(e) {
  const params = Object.fromEntries((e.params || []).map((p) => [p.key, p.value]));
  return e.poly ? e.poly() : e.zonohedron ? zonohedron(e.zonohedron(params)) : convexHull(e.solid(params));
}

function mount(card) {
  const e = CATALOG.find((x) => x.id === card.dataset.id);
  if (!e) return;
  const stage = card.querySelector('.mini-stage');
  const poly = polyOf(e);
  const inf = analyze(poly);
  const v = new Viewer(stage, { autoRotate: true, interactive: false });
  v.controls.autoRotateSpeed = 1.6;
  v.root.add(solidObject(poly, inf, { color: PALETTE[CATALOG.indexOf(e) % PALETTE.length] }));
  v.frame(inf.circumradius * 1.0);
}

export function setupGalleries() {
  const cards = [...document.querySelectorAll('.gallery .mini[data-id]')];
  for (const card of cards) {
    const e = CATALOG.find((x) => x.id === card.dataset.id);
    if (!e) continue;
    card.href = `catalog.html#${e.id}`;
    const note = card.dataset.note ? JSON.parse(card.dataset.note) : null;
    card.innerHTML = `<div class="mini-stage"></div><span class="mini-cap"><b>${L(e.name)}</b>${note ? `<span>${L(note)}</span>` : ''}</span>`;
  }
  const io = new IntersectionObserver((entries) => {
    for (const en of entries) {
      if (!en.isIntersecting) continue;
      io.unobserve(en.target);
      try { mount(en.target); } catch (err) { console.warn('mini3d', en.target.dataset.id, err); }
    }
  }, { rootMargin: '200px' });
  cards.forEach((c) => io.observe(c));
}
