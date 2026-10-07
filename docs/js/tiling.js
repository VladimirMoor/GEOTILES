// Построение фрагмента разбиения и его численная проверка.
import * as THREE from 'three';
import { v, det3, inverse3 } from './geom.js';

const mod2 = (x) => ((Math.round(x) % 2) + 2) % 2;

/** Параллелоэдр: решётка порождена векторами 2·(центр грани − центр тела). */
export function parallelohedronSpec(poly, info) {
  const vecs = poly.F.map((f) => v.mul(v.sub(v.avg(f.map((i) => poly.V[i])), info.center), 2));
  return { kind: 'translations', vecs };
}

function pickBasis(vecs) {
  let best = null, bestDet = Infinity;
  for (let i = 0; i < vecs.length; i++) for (let j = i + 1; j < vecs.length; j++) for (let k = j + 1; k < vecs.length; k++) {
    const d = Math.abs(det3(vecs[i], vecs[j], vecs[k]));
    if (d > 1e-9 && d < bestDet - 1e-9) { bestDet = d; best = [vecs[i], vecs[j], vecs[k]]; }
  }
  return best;
}

/**
 * spec — либо { kind:'translations', vecs }, либо { kind:'crystal', lattice:[a,b,c], motifs:[Matrix4] }.
 * Возвращает плитки с центрами не дальше R от центра исходной плитки.
 */
export function generateTiles(info, spec, R, maxTiles = 1500) {
  const c0 = new THREE.Vector3(...info.centroid);
  const tiles = [];
  let covolume;

  if (spec.kind === 'translations') {
    const basis = pickBasis(spec.vecs);
    covolume = Math.abs(det3(...basis));
    const coords = inverse3(basis);
    const seen = new Set(['0,0,0']);
    const queue = [[0, 0, 0]];
    while (queue.length && tiles.length < maxTiles) {
      const t = queue.shift();
      const [i, j, k] = coords(t);
      tiles.push({
        matrix: new THREE.Matrix4().makeTranslation(...t),
        centroid: c0.clone().add(new THREE.Vector3(...t)),
        cls: mod2(i) + 2 * mod2(j) + 4 * mod2(k),
        orbit: 0,
      });
      for (const f of spec.vecs) {
        const s = v.add(t, f);
        if (v.len(s) > R + 1e-9) continue;
        const key = s.map((x) => Math.round(x * 1e5)).join(',');
        if (!seen.has(key)) { seen.add(key); queue.push(s); }
      }
    }
  } else {
    const [a, b, c] = spec.lattice;
    covolume = Math.abs(det3(a, b, c)) / spec.motifs.length;
    const N = Math.min(10, Math.ceil(R / Math.min(v.len(a), v.len(b), v.len(c))) + 2);
    for (let i = -N; i <= N; i++) for (let j = -N; j <= N; j++) for (let k = -N; k <= N; k++) {
      const t = v.add(v.add(v.mul(a, i), v.mul(b, j)), v.mul(c, k));
      spec.motifs.forEach((m, mi) => {
        const M = new THREE.Matrix4().makeTranslation(...t).multiply(m);
        const cen = c0.clone().applyMatrix4(M);
        if (cen.distanceTo(c0) > R + 1e-9 || tiles.length >= maxTiles) return;
        tiles.push({ matrix: M, centroid: cen, cls: (mi * 3 + mod2(i) + 2 * mod2(j) + 4 * mod2(k)) % 8, orbit: mi });
      });
    }
  }
  return { tiles, covolume, c0 };
}

/**
 * Монте-Карло проверка: случайные точки внутри шара, который заведомо покрыт построенным фрагментом,
 * должны лежать строго внутри ровно одной плитки.
 */
export function checkCoverage(poly, info, gen, R, samples = 1500) {
  const r = R - info.circumradius;
  if (r <= 0.05) return null;
  const inv = gen.tiles.map((t) => t.matrix.clone().invert());
  const planes = poly.N.map((n, i) => ({ n, d: poly.D[i] }));
  const eps = 1e-9;
  let ok = 0, gaps = 0, overlaps = 0;
  const p = new THREE.Vector3(), q = new THREE.Vector3();
  for (let s = 0; s < samples; s++) {
    let x, y, z;
    do { x = Math.random() * 2 - 1; y = Math.random() * 2 - 1; z = Math.random() * 2 - 1; } while (x * x + y * y + z * z > 1);
    p.set(x * r, y * r, z * r).add(gen.c0);
    let hits = 0;
    for (let t = 0; t < inv.length && hits < 2; t++) {
      q.copy(p).applyMatrix4(inv[t]);
      const inside = planes.every(({ n, d }) => n[0] * q.x + n[1] * q.y + n[2] * q.z - d < -eps);
      if (inside) hits++;
    }
    if (hits === 1) ok++; else if (hits === 0) gaps++; else overlaps++;
  }
  return { samples, ok, gaps, overlaps, radius: r };
}
