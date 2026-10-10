// Три языка: английский (основной), русский, испанский.
// Язык: ?lang=xx → сохранённый выбор → en.

export const LANGS = [['en', 'EN'], ['ru', 'RU'], ['es', 'ES']];

function detect() {
  const q = new URLSearchParams(location.search).get('lang');
  if (LANGS.some(([c]) => c === q)) { try { localStorage.setItem('lang', q); } catch {} return q; }
  try { const s = localStorage.getItem('lang'); if (LANGS.some(([c]) => c === s)) return s; } catch {}
  return 'en';
}

export const lang = detect();
document.documentElement.lang = lang;

export function setLang(l) {
  try { localStorage.setItem('lang', l); } catch {}
  const u = new URL(location.href);
  u.searchParams.set('lang', l);
  location.href = u.toString();
}

/** Локализованное поле: строка или { en, ru, es }. */
export const L = (x) => (x == null || typeof x === 'string' ? x : x[lang] ?? x.en);

export function t(key, vars = {}) {
  const s = DICT[key] ? L(DICT[key]) : key;
  return s.replace(/\{(\w+)\}/g, (_, k) => vars[k] ?? '');
}

/** Подставляет переводы в элементы с data-i18n и показывает только блоки data-lang текущего языка. */
export function applyDOM(root = document) {
  root.querySelectorAll('[data-i18n]').forEach((el) => { el.innerHTML = t(el.dataset.i18n); });
  root.querySelectorAll('[data-lang]').forEach((el) => { el.hidden = el.dataset.lang !== lang; });
}

const NGONS = {
  3: { en: 'triangle', ru: 'треугольник', es: 'triángulo' },
  4: { en: 'quadrilateral', ru: 'четырёхугольник', es: 'cuadrilátero' },
  5: { en: 'pentagon', ru: 'пятиугольник', es: 'pentágono' },
  6: { en: 'hexagon', ru: 'шестиугольник', es: 'hexágono' },
  8: { en: 'octagon', ru: 'восьмиугольник', es: 'octágono' },
  10: { en: 'decagon', ru: 'десятиугольник', es: 'decágono' },
};

/** Название грани по коду из geom.classifyFace: 'square', 'rhombus', 'reg|5', 'equi|6', 'plain|7' … */
export function faceName(code) {
  const [kind, n] = code.split('|');
  if (!n) return t(`face.${kind}`);
  const name = NGONS[n] ? L(NGONS[n]) : t('face.ngon', { n });
  return kind === 'plain' ? name : t(`face.${kind}`, { name });
}

const DICT = {
  'nav.home': { en: 'Home', ru: 'Главная', es: 'Inicio' },
  'nav.catalog': { en: 'Catalog', ru: 'Каталог', es: 'Catálogo' },
  'nav.lab': { en: 'Lab', ru: 'Лаборатория', es: 'Laboratorio' },
  'nav.research': { en: 'Research', ru: 'Исследования', es: 'Investigación' },

  'title.home': { en: 'GEOTILES — open research on tilings', ru: 'GEOTILES — открытое исследование замощений', es: 'GEOTILES — investigación abierta sobre teselaciones' },
  'title.catalog': { en: 'Catalog · GEOTILES', ru: 'Каталог · GEOTILES', es: 'Catálogo · GEOTILES' },
  'title.lab': { en: 'Lab · GEOTILES', ru: 'Лаборатория · GEOTILES', es: 'Laboratorio · GEOTILES' },
  'title.research': { en: 'Research · GEOTILES', ru: 'Исследования · GEOTILES', es: 'Investigación · GEOTILES' },

  // главная
  'home.eyebrow': { en: 'open research · 2026', ru: 'открытое исследование · 2026', es: 'investigación abierta · 2026' },
  'home.h1': { en: 'Which shapes fill space — and how?', ru: 'Какие фигуры заполняют пространство — и как?', es: '¿Qué formas llenan el espacio, y cómo?' },
  'home.lead': {
    en: 'We work through a list of open tiling problems: equilateral polyhedra, tetrahedra, Heesch numbers, record space-fillers. Literature first, then exact computation. Everything is open: code, data, 3D models and preprints.',
    ru: 'Мы разбираем список открытых задач о замощениях: равносторонние многогранники, тетраэдры, числа Хееша, рекордные многогранники. Сначала литература, потом точные вычисления. Всё открыто: код, данные, 3D-модели и препринты.',
    es: 'Recorremos una lista de problemas abiertos de teselaciones: poliedros equiláteros, tetraedros, números de Heesch, poliedros récord. Primero la bibliografía, luego cálculo exacto. Todo es abierto: código, datos, modelos 3D y preprints.',
  },
  'home.cta.catalog': { en: 'Open the catalog', ru: 'Открыть каталог', es: 'Abrir el catálogo' },
  'home.cta.lab': { en: 'Lab', ru: 'Лаборатория', es: 'Laboratorio' },
  'home.cta.research': { en: 'Research', ru: 'Исследования', es: 'Investigación' },
  'home.tools.eyebrow': { en: 'tools', ru: 'инструменты', es: 'herramientas' },
  'home.tools.h': { en: 'What you can do here', ru: 'Что здесь можно делать', es: 'Qué puedes hacer aquí' },
  'home.card.catalog.h': { en: 'Catalog', ru: 'Каталог', es: 'Catálogo' },
  'home.card.catalog.p': {
    en: 'Equilateral space-fillers (parallelohedra, prisms, gyrobifastigia, our new 9- and 10-face solids) and record stereohedra, up to Engel’s 38-faced one. Each has a 3D model, a tiling fragment, a slice and an exploded view.',
    ru: 'Равносторонние тела, заполняющие пространство (параллелоэдры, призмы, гиробифастигиумы, наши новые 9- и 10-гранники), и рекордные стереоэдры вплоть до 38-гранника Энгеля. Для каждого — 3D-модель, фрагмент разбиения, срез и разлёт плиток.',
    es: 'Poliedros equiláteros que llenan el espacio (paraleloedros, prismas, girobifastigios, nuestros nuevos de 9 y 10 caras) y estereoedros récord, hasta el de 38 caras de Engel. Cada uno con modelo 3D, fragmento de teselación, corte y vista separada.',
  },
  'home.card.zono.h': { en: 'Zonohedron builder', ru: 'Конструктор зоноэдров', es: 'Constructor de zonoedros' },
  'home.card.zono.p': {
    en: 'Enter generators and get the solid. McMullen’s criterion tells whether it tiles space by translations, and the site builds the tiling.',
    ru: 'Задайте образующие и получите тело. Критерий МакМаллена скажет, замощает ли оно пространство сдвигами, а сайт построит это разбиение.',
    es: 'Introduce generadores y obtén el sólido. El criterio de McMullen indica si tesela el espacio por traslaciones, y el sitio construye la teselación.',
  },
  'home.card.hull.h': { en: 'Polyhedron analyzer', ru: 'Анализатор многогранников', es: 'Analizador de poliedros' },
  'home.card.hull.p': {
    en: 'Paste vertices. The site builds the hull, finds every dihedral angle and checks a necessary condition for tiling.',
    ru: 'Вставьте вершины. Сайт построит оболочку, найдёт все двугранные углы и проверит необходимое условие замощения.',
    es: 'Pega vértices. El sitio construye la envolvente, encuentra todos los ángulos diedros y comprueba una condición necesaria para teselar.',
  },
  'home.card.log.h': { en: 'Research', ru: 'Исследования', es: 'Investigación' },
  'home.card.log.p': {
    en: 'All problems with their status, results, preprints and notes, updated as the work goes on.',
    ru: 'Все задачи со статусом, результатами, препринтами и заметками; обновляется по ходу работы.',
    es: 'Todos los problemas con su estado, resultados, preprints y notas; se actualiza según avanza el trabajo.',
  },
  'home.problems.eyebrow': { en: 'problems', ru: 'задачи', es: 'problemas' },
  'home.problems.h': { en: 'Problem list', ru: 'Список задач', es: 'Lista de problemas' },
  'home.problems.p': {
    en: 'Ordered by our rough estimate of how realistic progress is, not by difficulty. Apart from the first, all are known open problems, so the realistic goal is modest: a new example, a new record or a closed special case.',
    ru: 'Задачи упорядочены по нашей оценке того, насколько реально продвинуться, а не по сложности. Кроме первой, все задачи известны и открыты, поэтому для них реалистичная цель скромнее: новый пример, новый рекорд или закрытый частный случай.',
    es: 'Ordenados según nuestra estimación de lo realista que es avanzar, no por dificultad. Salvo el primero, todos son problemas abiertos conocidos, así que el objetivo realista es modesto: un ejemplo nuevo, un récord nuevo o un caso particular resuelto.',
  },
  'footer': { en: 'GEOTILES · code and data are open on', ru: 'GEOTILES · код и данные открыты на', es: 'GEOTILES · código y datos abiertos en' },

  // панель инструментов
  'tb.solid': { en: 'Solid', ru: 'Тело', es: 'Sólido' },
  'tb.tiling': { en: 'Tiling', ru: 'Разбиение', es: 'Teselación' },
  'tb.mixed': { en: 'Mixed', ru: 'Пёстро', es: 'Mixto' },
  'tb.orbit': { en: 'By orientation', ru: 'По ориентации', es: 'Por orientación' },
  'tb.explode': { en: 'Explode', ru: 'Раздвинуть', es: 'Separar' },
  'tb.cut': { en: 'Slice', ru: 'Срез', es: 'Corte' },
  'tb.size': { en: 'Size', ru: 'Размер', es: 'Tamaño' },
  'tb.spin': { en: 'rotate', ru: 'вращать', es: 'girar' },

  'st.notile': { en: 'does not tile space', ru: 'не замощает пространство', es: 'no tesela el espacio' },
  'st.fragment': { en: 'tiling fragment', ru: 'фрагмент разбиения', es: 'fragmento de teselación' },
  'st.solid': { en: 'solid', ru: 'тело', es: 'sólido' },
  'st.faces': { en: '{n} faces', ru: 'граней: {n}', es: '{n} caras' },
  'st.equi': { en: 'equilateral', ru: 'равносторонний', es: 'equilátero' },
  'st.nonequi': { en: 'not equilateral', ru: 'не равносторонний', es: 'no equilátero' },

  'warn.outside': {
    en: 'For these parameters the union of the pieces is not convex: the convex hull is a different, non-equilateral solid. Move the sliders back into the family.',
    ru: 'При этих параметрах объединение частей невыпукло: выпуклая оболочка — другое, неравностороннее тело. Верните ползунки в пределы семейства.',
    es: 'Con estos parámetros la unión de las piezas no es convexa: la envolvente convexa es otro sólido, no equilátero. Vuelve a llevar los controles dentro de la familia.',
  },

  // свойства
  'p.params': { en: 'Parameters', ru: 'Параметры', es: 'Parámetros' },
  'p.comb': { en: 'Combinatorics', ru: 'Комбинаторика', es: 'Combinatoria' },
  'p.vef': { en: 'Vertices · edges · faces', ru: 'Вершины · рёбра · грани', es: 'Vértices · aristas · caras' },
  'p.euler': { en: 'Euler V − E + F', ru: 'Эйлер V − E + F', es: 'Euler V − E + F' },
  'p.sym': { en: 'Central symmetry', ru: 'Центральная симметрия', es: 'Simetría central' },
  'yes': { en: 'yes', ru: 'да', es: 'sí' },
  'no': { en: 'no', ru: 'нет', es: 'no' },
  'p.faces': { en: 'Faces', ru: 'Грани', es: 'Caras' },
  'p.metric': { en: 'Metric', ru: 'Метрика', es: 'Métrica' },
  'p.edges': { en: 'Edges, min … max', ru: 'Рёбра, min … max', es: 'Aristas, mín … máx' },
  'p.equi': { en: 'Equilateral', ru: 'Равносторонний', es: 'Equilátero' },
  'p.vol': { en: 'Volume', ru: 'Объём', es: 'Volumen' },
  'p.circ': { en: 'Circumradius', ru: 'Радиус описанной сферы', es: 'Radio circunscrito' },

  'face.tri': { en: 'triangle', ru: 'треугольник', es: 'triángulo' },
  'face.tri-reg': { en: 'equilateral triangle', ru: 'правильный треугольник', es: 'triángulo equilátero' },
  'face.square': { en: 'square', ru: 'квадрат', es: 'cuadrado' },
  'face.rhombus': { en: 'rhombus', ru: 'ромб', es: 'rombo' },
  'face.rect': { en: 'rectangle', ru: 'прямоугольник', es: 'rectángulo' },
  'face.quad': { en: 'quadrilateral', ru: 'четырёхугольник', es: 'cuadrilátero' },
  'face.reg': { en: 'regular {name}', ru: 'правильный {name}', es: '{name} regular' },
  'face.equi': { en: 'equilateral {name}', ru: 'равносторонний {name}', es: '{name} equilátero' },
  'face.ngon': { en: '{n}-gon', ru: '{n}-угольник', es: '{n}-ágono' },

  'f.title': { en: 'Dihedral angle filter', ru: 'Фильтр двугранных углов', es: 'Filtro de ángulos diedros' },
  'f.edges': { en: '{n} edges', ru: 'рёбер: {n}', es: '{n} aristas' },
  'f.fail': { en: 'cannot be completed to 360° or 180°', ru: 'не дополняется до 360° и 180°', es: 'no se completa a 360° ni a 180°' },
  'f.face': { en: '180° (neighbor’s face)', ru: '180° (грань соседа)', es: '180° (cara vecina)' },
  'f.ok': { en: 'Necessary condition holds: every angle is part of a 360° sum.', ru: 'Необходимое условие выполнено: каждый угол входит в сумму 360°.', es: 'La condición necesaria se cumple: cada ángulo forma parte de una suma de 360°.' },
  'f.bad': { en: 'Does not tile: around an edge with the marked angle, copies of this solid cannot close up 360°.', ru: 'Не замощает: вокруг ребра с отмеченным углом нельзя замкнуть 360° копиями этого тела.', es: 'No tesela: alrededor de una arista con el ángulo marcado, las copias de este sólido no pueden cerrar 360°.' },
  'f.help': {
    en: 'Around an interior point of an edge, neighboring tiles contribute either their dihedral angle or 180° (if the point lies inside their face), and the total is 360°.',
    ru: 'Вокруг внутренней точки ребра соседние плитки вносят свой двугранный угол или 180° (если точка лежит внутри их грани), и в сумме выходит 360°.',
    es: 'Alrededor de un punto interior de una arista, las teselas vecinas aportan su ángulo diedro o 180° (si el punto está dentro de su cara), y el total es 360°.',
  },

  'c.title': { en: 'Tiling check', ru: 'Проверка разбиения', es: 'Verificación de la teselación' },
  'c.tiles': { en: 'Tiles in fragment', ru: 'Плиток во фрагменте', es: 'Teselas en el fragmento' },
  'c.ratio': { en: 'V(tile) / V(fund. domain)', ru: 'V(плитки) / V(фунд. области)', es: 'V(tesela) / V(dominio fund.)' },
  'c.small': { en: 'Increase the fragment size to check coverage.', ru: 'Увеличьте размер фрагмента для проверки покрытия.', es: 'Aumenta el tamaño del fragmento para comprobar el recubrimiento.' },
  'c.cov': { en: '{ok} of {n} random points in a ball of radius {r} lie in exactly one tile', ru: '{ok} из {n} случайных точек в шаре радиуса {r} лежат ровно в одной плитке', es: '{ok} de {n} puntos aleatorios en una bola de radio {r} están en exactamente una tesela' },
  'c.covbad': { en: ' (gaps: {g}, overlaps: {o})', ru: ' (дыры: {g}, перекрытия: {o})', es: ' (huecos: {g}, solapamientos: {o})' },

  // лаборатория
  'lab.zono': { en: 'Zonohedron', ru: 'Зоноэдр', es: 'Zonoedro' },
  'lab.hull': { en: 'Polyhedron', ru: 'Многогранник', es: 'Poliedro' },
  'lab.zono.help': {
    en: 'Zonohedron generators, one per line: <code>x y z</code>. With normalization on, every edge has length 1 and the zonohedron is equilateral.',
    ru: 'Образующие зоноэдра, по одной на строку: <code>x y z</code>. Если включить нормировку, все рёбра получат длину 1 и зоноэдр станет равносторонним.',
    es: 'Generadores del zonoedro, uno por línea: <code>x y z</code>. Con la normalización activada, todas las aristas miden 1 y el zonoedro es equilátero.',
  },
  'lab.norm': { en: 'normalize generators', ru: 'нормировать образующие', es: 'normalizar generadores' },
  'lab.build': { en: 'Build', ru: 'Построить', es: 'Construir' },
  'lab.rand': { en: 'Random 4', ru: 'Случайные 4', es: '4 aleatorios' },
  'lab.hull.help': {
    en: 'Vertices, one per line: <code>x y z</code>. The convex hull is built and the dihedral angle filter is checked: if it fails, the solid certainly does not tile space.',
    ru: 'Вершины, по одной на строку: <code>x y z</code>. Строится выпуклая оболочка, затем проверяется фильтр двугранных углов: если он не выполнен, тело точно не замощает пространство.',
    es: 'Vértices, uno por línea: <code>x y z</code>. Se construye la envolvente convexa y se comprueba el filtro de ángulos diedros: si falla, el sólido seguro que no tesela el espacio.',
  },
  'lab.title.hull': { en: 'Convex hull', ru: 'Выпуклая оболочка', es: 'Envolvente convexa' },
  'm.title': { en: 'McMullen’s criterion', ru: 'Критерий МакМаллена', es: 'Criterio de McMullen' },
  'm.ok': { en: 'Parallelohedron: tiles space by translations', ru: 'Параллелоэдр: замощает пространство сдвигами', es: 'Paraleloedro: tesela el espacio por traslaciones' },
  'm.no': { en: 'Does not tile by translations', ru: 'Не замощает сдвигами', es: 'No tesela por traslaciones' },
  'm.r.ok': { en: 'along every generator at most 3 directions are seen', ru: 'вдоль любой образующей видно не больше 3 направлений', es: 'a lo largo de cada generador se ven como máximo 3 direcciones' },
  'm.r.planar': { en: 'generators are coplanar', ru: 'образующие лежат в одной плоскости', es: 'los generadores son coplanares' },
  'm.r.fail': { en: 'along generator ({g}) {k} directions are seen, more than 3', ru: 'вдоль образующей ({g}) видно {k} направлений, больше 3', es: 'a lo largo del generador ({g}) se ven {k} direcciones, más de 3' },
  'm.merged': { en: 'Generators after merging parallel ones: {n}.', ru: 'Образующих после слияния параллельных: {n}.', es: 'Generadores tras fusionar los paralelos: {n}.' },
  'm.rot': { en: 'A tiling with rotations is not excluded: see the angle filter below.', ru: 'Замощение с поворотами при этом не исключено: смотрите фильтр углов ниже.', es: 'No se descarta una teselación con rotaciones: mira el filtro de ángulos abajo.' },
  'err.points4': { en: 'At least 4 points are needed', ru: 'Нужно хотя бы 4 точки', es: 'Se necesitan al menos 4 puntos' },
  'err.coplanar': { en: 'The points are coplanar', ru: 'Точки лежат в одной плоскости', es: 'Los puntos son coplanares' },
  'err.gens8': { en: 'At most 8 generators', ru: 'Не больше 8 образующих', es: 'Como máximo 8 generadores' },
  'err.degenerate': { en: 'The solid is degenerate: some faces are almost coplanar', ru: 'Тело вырождено: некоторые грани почти лежат в одной плоскости', es: 'El sólido es degenerado: algunas caras son casi coplanares' },
  'err.parse': { en: 'Each line must contain three numbers', ru: 'Каждая строка должна содержать три числа', es: 'Cada línea debe contener tres números' },
};
