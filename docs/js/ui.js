// Общие куски интерфейса: шапка с переключателем языка и панели свойств многогранника.
import { t, lang, LANGS, setLang, applyDOM, faceName } from './i18n.js';

export const REPO = 'https://github.com/VladimirMoor/GEOTILES';

const LOGO = `<svg viewBox="0 0 32 32" aria-hidden="true"><g fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round">
<path d="M16 3 27 9.5v13L16 29 5 22.5v-13Z"/><path d="M16 3v13m0 0 11-6.5M16 16 5 9.5M16 16v13" opacity=".45"/></g></svg>`;

/** Шапка, заголовок страницы и переводы статических элементов. */
export function setupPage(active) {
  const links = [['index.html', 'home'], ['catalog.html', 'catalog'], ['lab.html', 'lab'], ['research.html', 'research']];
  const el = document.createElement('header');
  el.className = 'top';
  el.innerHTML = `<a class="brand" href="index.html">${LOGO}GEOTILES</a>
    <nav class="nav">${links.map(([href, id]) => `<a href="${href}" class="${id === active ? 'on' : ''}">${t(`nav.${id}`)}</a>`).join('')}</nav>
    <div class="langs">${LANGS.map(([c, label]) => `<button data-lang-btn="${c}" class="${c === lang ? 'on' : ''}">${label}</button>`).join('')}</div>
    <a class="gh" href="${REPO}" target="_blank" rel="noopener">GitHub ↗</a>`;
  document.body.prepend(el);
  el.querySelectorAll('[data-lang-btn]').forEach((b) => b.addEventListener('click', () => setLang(b.dataset.langBtn)));
  document.title = t(`title.${active}`);
  applyDOM();
}

const f = (x, d = 4) => Number(x).toFixed(d);
export const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

export function propsHTML(info) {
  const faces = Object.entries(info.faceTypes).map(([code, n]) => `<div>${n} × ${esc(faceName(code))}</div>`).join('');
  return `
    <h3>${t('p.comb')}</h3>
    <dl class="kv">
      <dt>${t('p.vef')}</dt><dd>${info.nV} · ${info.nE} · ${info.nF}</dd>
      <dt>${t('p.euler')}</dt><dd>${info.nV - info.nE + info.nF}</dd>
      <dt>${t('p.sym')}</dt><dd>${t(info.symmetric ? 'yes' : 'no')}</dd>
    </dl>
    <h3>${t('p.faces')}</h3><div class="small">${faces}</div>
    <h3>${t('p.metric')}</h3>
    <dl class="kv">
      <dt>${t('p.edges')}</dt><dd>${f(info.minLen, 6)} … ${f(info.maxLen, 6)}</dd>
      <dt>${t('p.equi')}</dt><dd><span class="chip ${info.equilateral ? 'ok' : 'bad'}">${t(info.equilateral ? 'yes' : 'no')}</span></dd>
      <dt>${t('p.vol')}</dt><dd>${f(info.volume, 6)}</dd>
      <dt>${t('p.circ')}</dt><dd>${f(info.circumradius)}</dd>
    </dl>`;
}

function exampleText(ex) {
  const parts = ex.terms.map(([k, a]) => (k > 1 ? `${k}×${a.toFixed(2)}°` : `${a.toFixed(2)}°`));
  if (ex.face) parts.push(t('f.face'));
  return `${parts.join(' + ')} = 360°`;
}

export function filterHTML(info) {
  const rows = info.filter.per.map((p) => {
    const cnt = info.dihedrals.find((d) => Math.abs(d.deg - p.deg) < 1e-6)?.count ?? '';
    return `<div class="angle ${p.ok ? '' : 'fail'}"><b>${p.deg.toFixed(2)}°</b>
      <span>${t('f.edges', { n: cnt })} · ${p.ok ? esc(exampleText(p.example)) : t('f.fail')}</span></div>`;
  }).join('');
  const verdict = `<div class="verdict ${info.filter.ok ? 'ok' : 'bad'}">${t(info.filter.ok ? 'f.ok' : 'f.bad')}</div>`;
  return `<h3>${t('f.title')}</h3><p class="small muted">${t('f.help')}</p><div class="angles">${rows}</div>${verdict}`;
}

export function tilingHTML(gen, info, cov) {
  if (!gen) return '';
  const ratio = info.volume / gen.covolume;
  const rOk = Math.abs(ratio - 1) < 1e-6;
  let covLine = `<div class="verdict dim">${t('c.small')}</div>`;
  if (cov) {
    const good = cov.ok === cov.samples;
    covLine = `<div class="verdict ${good ? 'ok' : 'bad'}">${good ? '✓' : '✗'} ${t('c.cov', { ok: cov.ok, n: cov.samples, r: cov.radius.toFixed(2) })}${good ? '' : t('c.covbad', { g: cov.gaps, o: cov.overlaps })}.</div>`;
  }
  return `<h3>${t('c.title')}</h3>
    <dl class="kv">
      <dt>${t('c.tiles')}</dt><dd>${gen.tiles.length}</dd>
      <dt>${t('c.ratio')}</dt><dd>${f(ratio, 6)} ${rOk ? '✓' : '✗'}</dd>
    </dl>${covLine}`;
}

export function mcmullenHTML(mc, merged) {
  const reason = mc.code === 'fail'
    ? t('m.r.fail', { g: mc.g.map((x) => x.toFixed(3)).join(', '), k: mc.k })
    : t(`m.r.${mc.code}`);
  return `<h3>${t('m.title')}</h3>
    <div class="verdict ${mc.ok ? 'ok' : 'bad'}">${t(mc.ok ? 'm.ok' : 'm.no')}: ${reason}.</div>
    <p class="small muted" style="margin-top:8px">${t('m.merged', { n: merged })} ${mc.ok ? '' : t('m.rot')}</p>`;
}

/** Полоска вкладок задач на страницах исследований: «Все задачи · 01 · 02 · 03 · 04». */
export async function problemTabs(active) {
  const el = document.getElementById('ptabs');
  if (!el) return;
  const { PROBLEMS, PAGES } = await import('./problems.js');
  const { L } = await import('./i18n.js');
  const all = { en: 'All problems', ru: 'Все задачи', es: 'Todos los problemas' };
  if (active > 0) document.title = `${String(active).padStart(2, '0')} · ${L(PROBLEMS[active - 1][0])} · GEOTILES`;
  el.innerHTML = `<a href="research.html" class="${active === 0 ? 'on' : ''}">${L(all)}</a>` +
    PAGES.map((href, i) => `<a href="${href}" class="${active === i + 1 ? 'on' : ''}"><b>${String(i + 1).padStart(2, '0')}</b> ${L(PROBLEMS[i][0])}</a>`).join('');
}
