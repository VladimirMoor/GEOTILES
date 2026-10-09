import { zonohedron, analyze } from './geom.js';
import { Viewer, tilingObject } from './viewer.js';
import { generateTiles, parallelohedronSpec } from './tiling.js';
import { setupPage } from './ui.js';
import { L } from './i18n.js';

setupPage('home');

const S = {
  work: { en: 'in progress', ru: 'в работе', es: 'en curso' },
  open: { en: 'open', ru: 'открыта', es: 'abierto' },
  far: { en: 'out of reach now', ru: 'пока недостижимо', es: 'fuera de alcance' },
  fund: { en: 'fundamental', ru: 'фундаментальная', es: 'fundamental' },
};

const PROBLEMS = [
  [{ en: 'Equilateral convex polyhedra that tile space', ru: 'Равносторонние выпуклые многогранники, замощающие пространство', es: 'Poliedros convexos equiláteros que teselan el espacio' },
    { en: 'Our own question; no classification exists in the literature. First pass for ≤ 8 faces is done: all 301 combinatorial types processed, no rigid solid tiles, 11 families excluded by the Dehn invariant, 5 families tile in full. Pentagonal and hexagonal prisms remain open.', ru: 'Наш собственный вопрос; в литературе классификации нет. Первый проход для ≤ 8 граней готов: обработаны все 301 комбинаторный тип, ни одно жёсткое тело не замощает, 11 семейств исключены инвариантом Дена, 5 семейств замощают целиком. Пятиугольные и шестиугольные призмы пока открыты.', es: 'Nuestra propia pregunta; no existe clasificación en la bibliografía. La primera pasada para ≤ 8 caras está hecha: procesados los 301 tipos combinatorios, ningún sólido rígido tesela, 11 familias descartadas por el invariante de Dehn y 5 familias teselan por completo. Los prismas pentagonales y hexagonales siguen abiertos.' }, 'work', 'info'],
  [{ en: 'New space-filling tetrahedra', ru: 'Новые тетраэдры, замощающие пространство', es: 'Nuevos tetraedros que llenan el espacio' },
    { en: 'Target: the 40 sporadic rational tetrahedra left undecided by Chentouf–Sun (2023). Two of them turn out to be Sommerville tilers; none of the other 38 tiles face-to-face. Non-face-to-face tilings remain open.', ru: 'Цель — 40 спорадических рациональных тетраэдров, оставленных нерешёнными у Chentouf–Sun (2023). Два из них оказались замостителями Соммервилля; ни один из остальных 38 не замощает «лицом к лицу». Разбиения не «лицом к лицу» пока открыты.', es: 'Objetivo: los 40 tetraedros racionales esporádicos que Chentouf–Sun (2023) dejaron sin decidir. Dos resultan ser teselas de Sommerville; ninguno de los otros 38 tesela cara a cara. Las teselaciones no cara a cara siguen abiertas.' }, 'work', 'info'],
  [{ en: 'A polygon with Heesch number 7', ru: 'Многоугольник с числом Хееша 7', es: 'Un polígono con número de Heesch 7' },
    { en: 'Record is 6 (Bašić 2021). We built Kaplan’s SAT tool and ran the first Heesch census of polydrafters: all 9.7 million with ≤ 20 cells; the maximum is 3. Next: families around the known Heesch-6 tiles.', ru: 'Рекорд — 6 (Башич, 2021). Мы собрали SAT-инструмент Каплана и провели первую перепись чисел Хееша для полидрафтеров: все 9,7 млн фигур до 20 клеток; максимум — 3. Дальше — семейства вокруг известных плиток с числом 6.', es: 'El récord es 6 (Bašić, 2021). Compilamos la herramienta SAT de Kaplan e hicimos el primer censo de números de Heesch de polidrafters: los 9,7 millones con ≤ 20 celdas; el máximo es 3. Siguiente: familias alrededor de las teselas conocidas con número 6.' }, 'work', 'info'],
  [{ en: 'A space-filler with more than 38 faces', ru: 'Больше 38 граней у многогранника, замощающего пространство', es: 'Un poliedro que llena el espacio con más de 38 caras' },
    { en: 'Another record: search among stereohedra of crystallographic groups.', ru: 'Тоже рекорд: поиск среди стереоэдров кристаллографических групп.', es: 'Otro récord: búsqueda entre estereoedros de grupos cristalográficos.' }, 'open', 'dim'],
  [{ en: 'Periodic tiling conjecture in dimension 3', ru: 'Периодичность в размерности 3', es: 'Conjetura de teselación periódica en dimensión 3' },
    { en: 'Must a tile that tiles space by translations also admit a periodic tiling? (Greenfeld–Tao)', ru: 'Обязательно ли плитка, замощающая пространство сдвигами, допускает и периодическое замощение? (Гринфельд–Тао)', es: '¿Una tesela que llena el espacio por traslaciones admite siempre una teselación periódica? (Greenfeld–Tao)' }, 'open', 'dim'],
  [{ en: 'Voronoi’s conjecture in dimension ≥ 6', ru: 'Гипотеза Вороного для размерности ≥ 6', es: 'Conjetura de Voronói en dimensión ≥ 6' },
    { en: 'Needs deep theory; enumeration will not do it.', ru: 'Нужна глубокая теория, перебором её не взять.', es: 'Requiere teoría profunda; la enumeración no basta.' }, 'open', 'dim'],
  [{ en: 'Classifying convex polyhedra that tile space', ru: 'Классификация выпуклых многогранников, замощающих пространство', es: 'Clasificar los poliedros convexos que teselan el espacio' },
    { en: 'The 3D analogue of Rao’s result.', ru: 'Трёхмерный аналог результата Рао.', es: 'El análogo 3D del resultado de Rao.' }, 'far', 'dim'],
  [{ en: 'Decidability of tiling by one polygon', ru: 'Разрешимость замощения одним многоугольником', es: 'Decidibilidad de teselar con un polígono' },
    { en: 'Is there an algorithm that decides whether a given polygon tiles the plane?', ru: 'Есть ли алгоритм, который по многоугольнику решает, замощает ли он плоскость?', es: '¿Existe un algoritmo que decida si un polígono dado tesela el plano?' }, 'fund', 'dim'],
  [{ en: 'Are Heesch numbers bounded?', ru: 'Ограничены ли числа Хееша?', es: '¿Están acotados los números de Heesch?' },
    { en: 'The simplest statement on the list and perhaps the hardest question.', ru: 'Самая простая формулировка и, возможно, самый трудный вопрос в списке.', es: 'El enunciado más simple de la lista y quizá la pregunta más difícil.' }, 'fund', 'dim'],
];

const LINKS = ['research.html', 'tetrahedra.html'];
document.getElementById('problems').innerHTML = PROBLEMS.map(([title, desc, s, c], i) => `
  <div class="problem ${LINKS[i] ? 'active' : ''}"><div class="n">${String(i + 1).padStart(2, '0')}</div>
  <div><h3>${LINKS[i] ? `<a href="${LINKS[i]}">${L(title)}</a>` : L(title)}</h3><p>${L(desc)}</p></div><span class="chip ${c}">${L(S[s])}</span></div>`).join('');

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
