// Вычислительная геометрия для выпуклых многогранников: оболочка, рёбра, двугранные углы,
// зоноэдры, критерий МакМаллена и необходимый фильтр двугранных углов.

export const v = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
  mul: (a, s) => [a[0] * s, a[1] * s, a[2] * s],
  dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
  cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
  len: (a) => Math.hypot(a[0], a[1], a[2]),
  unit: (a) => { const l = Math.hypot(a[0], a[1], a[2]); return [a[0] / l, a[1] / l, a[2] / l]; },
  avg: (pts) => v.mul(pts.reduce((s, p) => v.add(s, p), [0, 0, 0]), 1 / pts.length),
};

export const det3 = (a, b, c) => v.dot(a, v.cross(b, c));
export const DEG = 180 / Math.PI;

export function inverse3(cols) {
  // cols — три вектора-столбца; возвращает функцию p -> координаты p в этом базисе
  const [a, b, c] = cols;
  const d = det3(a, b, c);
  const r0 = v.mul(v.cross(b, c), 1 / d), r1 = v.mul(v.cross(c, a), 1 / d), r2 = v.mul(v.cross(a, b), 1 / d);
  return (p) => [v.dot(r0, p), v.dot(r1, p), v.dot(r2, p)];
}

/** Выпуклая оболочка небольшого набора точек (перебор троек, O(n^4)). Грани — многоугольники, обход против часовой стрелки снаружи. */
export function convexHull(input) {
  const scale = Math.max(1e-9, ...input.map((p) => v.len(p)));
  const tol = 1e-7 * scale;
  const pts = [];
  for (const p of input) if (!pts.some((q) => v.len(v.sub(p, q)) < tol * 10)) pts.push(p);
  const n = pts.length;
  if (n < 4) throw new Error('err.points4');

  const planes = [];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) for (let k = j + 1; k < n; k++) {
    let nr = v.cross(v.sub(pts[j], pts[i]), v.sub(pts[k], pts[i]));
    const L = v.len(nr);
    if (L < tol) continue;
    nr = v.mul(nr, 1 / L);
    let d = v.dot(nr, pts[i]);
    let pos = 0, neg = 0;
    for (const p of pts) {
      const s = v.dot(nr, p) - d;
      if (s > tol) pos++; else if (s < -tol) neg++;
      if (pos && neg) break;
    }
    if ((pos && neg) || (!pos && !neg)) continue;
    if (pos) { nr = v.mul(nr, -1); d = -d; }
    if (!planes.some((q) => v.dot(q.n, nr) > 1 - 1e-9 && Math.abs(q.d - d) < tol * 10)) planes.push({ n: nr, d });
  }
  if (planes.length < 4) throw new Error('err.coplanar');

  const faces = [];
  for (const { n: nr, d } of planes) {
    const onPlane = [];
    pts.forEach((p, i) => { if (Math.abs(v.dot(nr, p) - d) < tol * 10) onPlane.push(i); });
    // плоская выпуклая оболочка точек грани: отбрасывает точки внутри грани и на рёбрах
    const c = v.avg(onPlane.map((i) => pts[i]));
    const u = v.unit(v.sub(pts[onPlane[0]], c));
    const w = v.cross(nr, u);
    const q = onPlane.map((i) => { const r = v.sub(pts[i], c); return { i, x: v.dot(r, u), y: v.dot(r, w) }; })
      .sort((a, b) => a.x - b.x || a.y - b.y);
    const turn = (o, a, b) => (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    const chain = (list) => {
      const out = [];
      for (const p of list) {
        while (out.length >= 2 && turn(out[out.length - 2], out[out.length - 1], p) <= tol * scale) out.pop();
        out.push(p);
      }
      out.pop();
      return out;
    };
    const idx = [...chain(q), ...chain([...q].reverse())].map((p) => p.i);
    faces.push({ idx, n: nr, d });
  }

  const used = [...new Set(faces.flatMap((f) => f.idx))];
  const remap = new Map(used.map((old, i) => [old, i]));
  return {
    V: used.map((i) => pts[i]),
    F: faces.map((f) => f.idx.map((i) => remap.get(i))),
    N: faces.map((f) => f.n),
    D: faces.map((f) => f.d),
  };
}

/** Код типа грани для i18n.faceName: 'tri-reg', 'square', 'rhombus', 'reg|5', 'equi|6', 'plain|7' … */
function classifyFace(pts) {
  const k = pts.length;
  const lens = pts.map((p, i) => v.len(v.sub(pts[(i + 1) % k], p)));
  const equi = Math.max(...lens) / Math.min(...lens) - 1 < 1e-6;
  const angs = pts.map((p, i) => {
    const a = v.sub(pts[(i - 1 + k) % k], p), b = v.sub(pts[(i + 1) % k], p);
    return Math.acos(Math.max(-1, Math.min(1, v.dot(a, b) / (v.len(a) * v.len(b)))));
  });
  const equiAng = Math.max(...angs) - Math.min(...angs) < 1e-6;
  if (k === 3) return equi ? 'tri-reg' : 'tri';
  if (k === 4) {
    if (equi) return equiAng ? 'square' : 'rhombus';
    return equiAng ? 'rect' : 'quad';
  }
  if (equi && equiAng) return `reg|${k}`;
  return equi ? `equi|${k}` : `plain|${k}`;
}

/** Все свойства, которые показываются в карточке многогранника. */
export function analyze(poly) {
  const { V, F, N } = poly;
  const edgeMap = new Map();
  F.forEach((f, fi) => f.forEach((a, t) => {
    const b = f[(t + 1) % f.length];
    const key = a < b ? `${a}-${b}` : `${b}-${a}`;
    if (!edgeMap.has(key)) edgeMap.set(key, { a: Math.min(a, b), b: Math.max(a, b), faces: [] });
    edgeMap.get(key).faces.push(fi);
  }));
  const edges = [...edgeMap.values()].map((e) => {
    const [f1, f2] = e.faces;
    const cos = Math.max(-1, Math.min(1, v.dot(N[f1], N[f2])));
    return { ...e, length: v.len(v.sub(V[e.a], V[e.b])), dihedral: Math.PI - Math.acos(cos) };
  });

  const lens = edges.map((e) => e.length);
  const minLen = Math.min(...lens), maxLen = Math.max(...lens);

  // объём и центр масс разбиением на тетраэдры из внутренней точки
  const o = v.avg(V);
  let vol = 0, cm = [0, 0, 0];
  for (const f of F) for (let i = 1; i < f.length - 1; i++) {
    const a = v.sub(V[f[0]], o), b = v.sub(V[f[i]], o), c = v.sub(V[f[i + 1]], o);
    const tv = det3(a, b, c) / 6;
    vol += tv;
    cm = v.add(cm, v.mul(v.add(v.add(a, b), c), tv / 4));
  }
  const centroid = v.add(o, v.mul(cm, 1 / vol));
  const circumradius = Math.max(...V.map((p) => v.len(v.sub(p, centroid))));

  const faceTypes = {};
  for (const f of F) {
    const t = classifyFace(f.map((i) => V[i]));
    faceTypes[t] = (faceTypes[t] || 0) + 1;
  }

  const dihedrals = [];
  for (const e of edges) {
    const deg = e.dihedral * DEG;
    const hit = dihedrals.find((d) => Math.abs(d.deg - deg) < 1e-6);
    if (hit) hit.count++; else dihedrals.push({ deg, count: 1 });
  }
  dihedrals.sort((a, b) => a.deg - b.deg);

  const symmetric = V.every((p) => {
    const q = v.sub(v.mul(o, 2), p);
    return V.some((r) => v.len(v.sub(r, q)) < 1e-6 * circumradius);
  });

  return {
    nV: V.length, nE: edges.length, nF: F.length,
    edges, minLen, maxLen,
    equilateral: maxLen / minLen - 1 < 1e-6,
    volume: vol, centroid, center: o, circumradius,
    faceTypes, dihedrals, symmetric,
    filter: dihedralFilter(dihedrals.map((d) => d.deg)),
  };
}

/**
 * Необходимое условие монотайлинга: вокруг внутренней точки каждого ребра плитки углы соседних плиток
 * дают 360°. Каждая соседняя плитка вносит свой двугранный угол или 180° (если точка лежит внутри её грани),
 * и 180° может быть не больше одного раза. Значит, каждый угол α должен входить в сумму
 * Σ kᵢ·αᵢ = 360° или Σ kᵢ·αᵢ = 180° с неотрицательными целыми kᵢ.
 */
export function dihedralFilter(angles, tol = 1e-5) {
  const A = [...angles].sort((a, b) => b - a);
  const solve = (target) => {
    const sols = [];
    const k = new Array(A.length).fill(0);
    const rec = (i, rest) => {
      if (sols.length >= 400) return;
      if (Math.abs(rest) < tol) { sols.push([...k]); return; }
      if (i >= A.length || rest < -tol) return;
      for (let m = Math.floor((rest + tol) / A[i]); m >= 0; m--) {
        k[i] = m;
        rec(i + 1, rest - m * A[i]);
      }
      k[i] = 0;
    };
    rec(0, target);
    return sols.filter((s) => s.some((x) => x > 0));
  };
  const s360 = solve(360), s180 = solve(180);
  const per = A.map((deg, i) => {
    const full = s360.find((s) => s[i] > 0);
    const half = s180.find((s) => s[i] > 0);
    const terms = (s) => A.flatMap((a, j) => (s[j] ? [[s[j], a]] : []));
    return {
      deg, ok: !!(full || half),
      example: full ? { terms: terms(full), face: false } : half ? { terms: terms(half), face: true } : null,
    };
  }).sort((a, b) => a.deg - b.deg);
  return { ok: per.every((p) => p.ok), per };
}

/** Убирает нулевые и сливает параллельные образующие. */
export function mergeGenerators(gens) {
  const out = [];
  for (const g of gens) {
    if (v.len(g) < 1e-9) continue;
    const u = v.unit(g);
    const j = out.findIndex((h) => Math.abs(Math.abs(v.dot(v.unit(h), u)) - 1) < 1e-9);
    if (j < 0) out.push(g);
    else out[j] = v.add(out[j], v.dot(out[j], g) > 0 ? g : v.mul(g, -1));
  }
  return out;
}

/** Зоноэдр Σ[−gᵢ/2, gᵢ/2]. */
export function zonohedron(gens) {
  const G = mergeGenerators(gens);
  if (G.length > 8) throw new Error('err.gens8');
  const pts = [];
  for (let m = 0; m < 1 << G.length; m++) {
    let p = [0, 0, 0];
    G.forEach((g, i) => { p = v.add(p, v.mul(g, (m >> i) & 1 ? 0.5 : -0.5)); });
    pts.push(p);
  }
  return convexHull(pts);
}

/**
 * Критерий МакМаллена (1975) в размерности 3: зоноэдр замощает пространство параллельными переносами тогда
 * и только тогда, когда для каждой образующей g проекции остальных на плоскость ⊥ g дают не больше трёх направлений.
 */
export function mcmullen(gens) {
  const G = mergeGenerators(gens);
  const rank3 = G.some((a, i) => G.some((b, j) => j > i && G.some((c, k) => k > j && Math.abs(det3(v.unit(a), v.unit(b), v.unit(c))) > 1e-9)));
  if (!rank3) return { ok: false, code: 'planar' };
  let worst = 0, worstG = null;
  for (const g of G) {
    const u = v.unit(g);
    const dirs = [];
    for (const h of G) {
      const p = v.sub(h, v.mul(u, v.dot(h, u)));
      if (v.len(p) < 1e-9) continue;
      const pu = v.unit(p);
      if (!dirs.some((d) => Math.abs(Math.abs(v.dot(d, pu)) - 1) < 1e-9)) dirs.push(pu);
    }
    if (dirs.length > worst) { worst = dirs.length; worstG = g; }
  }
  return worst <= 3 ? { ok: true, code: 'ok' } : { ok: false, code: 'fail', k: worst, g: worstG };
}
