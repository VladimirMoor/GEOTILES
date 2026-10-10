// Каталог равносторонних выпуклых многогранников (длина ребра 1). Тексты на трёх языках.
// zonohedron(p) → образующие, или solid(p) → вершины; tiling(p) → 'parallelohedron' | { lattice, motifs } | нет.
import * as THREE from 'three';
import { STEREO } from './stereo-data.js';

const h = Math.sqrt(3) / 2;
const phi = (1 + Math.sqrt(5)) / 2;
const rad = (d) => (d * Math.PI) / 180;
const unit = (a) => { const l = Math.hypot(...a); return a.map((x) => x / l); };
const prism = (base, H = 1) => [...base.map(([x, y]) => [x, y, 0]), ...base.map(([x, y]) => [x, y, H])];
const rotZ = (a, t = [0, 0, 0]) => new THREE.Matrix4().makeRotationZ(a).setPosition(...t);
const I = () => new THREE.Matrix4();

const opMatrix = ([R, t]) => new THREE.Matrix4().set(R[0][0], R[0][1], R[0][2], t[0], R[1][0], R[1][1], R[1][2], t[1], R[2][0], R[2][1], R[2][2], t[2], 0, 0, 0, 1);
const stereo = (key) => ({ poly: () => ({ V: STEREO[key].points, F: STEREO[key].faces, N: STEREO[key].normals, D: STEREO[key].offsets }), solid: () => STEREO[key].points, tiling: () => ({ lattice: STEREO[key].lattice, motifs: STEREO[key].ops.map(opMatrix) }) });

export const GROUPS = [
  { id: 'par', title: { en: 'Parallelohedra', ru: 'Параллелоэдры', es: 'Paraleloedros' },
    note: { en: 'Tile space by translations alone', ru: 'Замощают пространство параллельными переносами', es: 'Teselan el espacio solo por traslaciones' } },
  { id: 'prism', title: { en: 'Prisms', ru: 'Призмы', es: 'Prismas' },
    note: { en: 'The base tiles the plane; layers stack up', ru: 'Основание замощает плоскость, слои укладываются стопкой', es: 'La base tesela el plano y las capas se apilan' } },
  { id: 'stereo', title: { en: 'Other space-fillers', ru: 'Другие замощающие', es: 'Otros que llenan el espacio' },
    note: { en: 'Rotations are needed', ru: 'Нужны повороты', es: 'Se necesitan rotaciones' } },
  { id: 'found', title: { en: 'Found by our search', ru: 'Найдены нашим перебором', es: 'Encontrados por nuestra búsqueda' },
    note: { en: 'Continuous families of equilateral space-fillers with ≤ 8 faces; every member tiles', ru: 'Непрерывные семейства равносторонних тел с ≤ 8 гранями; замощает каждый член', es: 'Familias continuas de poliedros equiláteros con ≤ 8 caras; todos sus miembros teselan' } },
  { id: 'records', title: { en: 'Record stereohedra', ru: 'Рекордные стереоэдры', es: 'Estereoedros récord' },
    note: { en: 'Voronoi cells of symmetric point sets with the most faces; not equilateral. Facet counts certified in exact arithmetic', ru: 'Ячейки Вороного симметричных наборов точек с наибольшим числом граней; не равносторонние. Число граней подтверждено точной арифметикой', es: 'Celdas de Voronói de conjuntos simétricos de puntos con más caras; no equiláteras. Número de caras certificado con aritmética exacta' } },
  { id: 'no', title: { en: 'Do not tile', ru: 'Не замощают', es: 'No teselan' },
    note: { en: 'For comparison: ruled out by the dihedral angle filter', ru: 'Для сравнения: исключены фильтром двугранных углов', es: 'Para comparar: descartados por el filtro de ángulos diedros' } },
];

export const CATALOG = [
  {
    id: 'cube', group: 'par', status: 'tiles',
    name: { en: 'Cube', ru: 'Куб', es: 'Cubo' },
    zonohedron: () => [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tiling: () => 'parallelohedron',
    about: {
      en: 'The only Platonic solid that fills space.',
      ru: 'Единственное платоново тело, заполняющее пространство.',
      es: 'El único sólido platónico que llena el espacio.',
    },
  },
  {
    id: 'rhombohedron', group: 'par', status: 'tiles',
    name: { en: 'Rhombohedron', ru: 'Ромбоэдр', es: 'Romboedro' },
    params: [{ key: 'theta', label: { en: 'angle between edges', ru: 'угол между рёбрами', es: 'ángulo entre aristas' }, min: 20, max: 118, step: 0.5, value: 70, unit: '°' }],
    zonohedron: ({ theta }) => {
      const c = Math.cos(rad(theta)), s = Math.sin(rad(theta));
      const cy = (c - c * c) / s;
      return [[1, 0, 0], [c, s, 0], [c, cy, Math.sqrt(Math.max(0, 1 - c * c - cy * cy))]];
    },
    tiling: () => 'parallelohedron',
    about: {
      en: 'A parallelepiped whose faces are all congruent rhombi. Every parallelepiped tiles space by translations, so the whole family does. At 90° it is the cube. A general equilateral parallelepiped (three arbitrary unit generators) can be built in the Lab.',
      ru: 'Параллелепипед, у которого все грани — равные ромбы. Любой параллелепипед замощает пространство сдвигами, поэтому это семейство целиком. При 90° получается куб. Общий равносторонний параллелепипед (три произвольных единичных образующих) соберите в лаборатории.',
      es: 'Un paralelepípedo cuyas caras son rombos congruentes. Todo paralelepípedo tesela el espacio por traslaciones, así que toda la familia lo hace. Con 90° es el cubo. Un paralelepípedo equilátero general (tres generadores unitarios arbitrarios) se puede construir en el Laboratorio.',
    },
  },
  {
    id: 'hex-prism', group: 'par', status: 'tiles',
    name: { en: 'Hexagonal prism', ru: 'Шестиугольная призма', es: 'Prisma hexagonal' },
    params: [
      { key: 'beta', label: { en: 'direction of 2nd edge', ru: 'направление 2-го ребра', es: 'dirección de la 2.ª arista' }, min: 15, max: 150, step: 0.5, value: 60, unit: '°' },
      { key: 'gamma', label: { en: 'direction of 3rd edge', ru: 'направление 3-го ребра', es: 'dirección de la 3.ª arista' }, min: 30, max: 165, step: 0.5, value: 120, unit: '°' },
    ],
    zonohedron: ({ beta, gamma }) => {
      const b = Math.min(beta, gamma - 5), g = Math.max(gamma, beta + 5);
      return [[1, 0, 0], [Math.cos(rad(b)), Math.sin(rad(b)), 0], [Math.cos(rad(g)), Math.sin(rad(g)), 0], [0, 0, 1]];
    },
    tiling: () => 'parallelohedron',
    about: {
      en: 'A prism over a centrally symmetric equilateral hexagon with square sides. At 60°/120° it is the regular hexagonal prism. A two-parameter family.',
      ru: 'Призма над центрально-симметричным равносторонним шестиугольником с квадратными боковыми гранями. При 60°/120° — правильная шестиугольная призма. Двухпараметрическое семейство.',
      es: 'Un prisma sobre un hexágono equilátero centralmente simétrico con caras laterales cuadradas. Con 60°/120° es el prisma hexagonal regular. Una familia de dos parámetros.',
    },
  },
  {
    id: 'rhombic-dodecahedron', group: 'par', status: 'tiles',
    name: { en: 'Kepler’s rhombic dodecahedron', ru: 'Ромбододекаэдр Кеплера', es: 'Dodecaedro rómbico de Kepler' },
    zonohedron: () => [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(unit),
    tiling: () => 'parallelohedron',
    about: {
      en: '12 congruent rhombi with diagonal ratio √2: the Voronoi cell of the face-centred cubic lattice. Any zonohedron with four unit generators in general position is an equilateral “rhombic dodecahedron”, and all of them tile by translations — a five-parameter family.',
      ru: '12 равных ромбов с отношением диагоналей √2. Ячейка Вороного гранецентрированной кубической решётки. Любой зоноэдр из четырёх единичных образующих общего положения — равносторонний «ромбододекаэдр», и все они замощают сдвигами: это пятипараметрическое семейство.',
      es: '12 rombos congruentes con razón de diagonales √2: la celda de Voronoi de la red cúbica centrada en las caras. Todo zonoedro con cuatro generadores unitarios en posición general es un «dodecaedro rómbico» equilátero, y todos teselan por traslaciones: una familia de cinco parámetros.',
    },
  },
  {
    id: 'bilinski', group: 'par', status: 'tiles',
    name: { en: 'Bilinski dodecahedron', ru: 'Додекаэдр Билинского', es: 'Dodecaedro de Bilinski' },
    zonohedron: () => [[0, 1, phi], [0, -1, phi], [1, phi, 0], [-1, phi, 0]].map(unit),
    tiling: () => 'parallelohedron',
    about: {
      en: '12 golden rhombi, Bilinski (1960). Four of the six fivefold axes of the icosahedron. Combinatorially a rhombic dodecahedron, and also a parallelohedron.',
      ru: '12 золотых ромбов, Билинский (1960). Четыре из шести осей пятого порядка икосаэдра. Комбинаторно эквивалентен ромбододекаэдру и тоже параллелоэдр.',
      es: '12 rombos áureos, Bilinski (1960). Cuatro de los seis ejes de orden cinco del icosaedro. Combinatoriamente es un dodecaedro rómbico y también un paraleloedro.',
    },
  },
  {
    id: 'elongated-dodecahedron', group: 'par', status: 'tiles',
    name: { en: 'Elongated dodecahedron', ru: 'Удлинённый додекаэдр', es: 'Dodecaedro alargado' },
    zonohedron: () => [...[[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(unit), [0, 0, 1]],
    tiling: () => 'parallelohedron',
    about: {
      en: '8 rhombi and 4 equilateral hexagons. The fourth of Fedorov’s parallelohedron types.',
      ru: '8 ромбов и 4 равносторонних шестиугольника. Четвёртый тип параллелоэдров Фёдорова.',
      es: '8 rombos y 4 hexágonos equiláteros. El cuarto tipo de paraleloedro de Fiódorov.',
    },
  },
  {
    id: 'truncated-octahedron', group: 'par', status: 'tiles',
    name: { en: 'Truncated octahedron', ru: 'Усечённый октаэдр', es: 'Octaedro truncado' },
    zonohedron: () => [[1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]].map(unit),
    tiling: () => 'parallelohedron',
    about: {
      en: '6 squares and 8 regular hexagons: the Voronoi cell of the body-centred lattice, Kelvin’s foam. The fifth and most complex Fedorov type.',
      ru: '6 квадратов и 8 правильных шестиугольников. Ячейка Вороного объёмно-центрированной решётки, «пена Кельвина». Пятый тип Фёдорова, самый сложный.',
      es: '6 cuadrados y 8 hexágonos regulares: la celda de Voronoi de la red centrada en el cuerpo, la espuma de Kelvin. El quinto tipo de Fiódorov, el más complejo.',
    },
  },
  {
    id: 'tri-prism', group: 'prism', status: 'tiles',
    name: { en: 'Triangular prism', ru: 'Треугольная призма', es: 'Prisma triangular' },
    solid: () => prism([[0, 0], [1, 0], [0.5, h]]),
    tiling: () => ({ lattice: [[1, 0, 0], [0.5, h, 0], [0, 0, 1]], motifs: [I(), rotZ(Math.PI, [1.5, h, 0])] }),
    about: {
      en: 'Regular-faced. Equilateral triangular prisms are exactly the prisms with an equilateral triangle base and rhombic sides, and all of them tile.',
      ru: 'Правильногранная. Равносторонние треугольные призмы — это все призмы с правильным треугольником в основании и ромбическими боковыми гранями; все они замощают.',
      es: 'De caras regulares. Los prismas triangulares equiláteros son exactamente los de base triángulo equilátero y caras laterales rómbicas, y todos teselan.',
    },
  },
  {
    id: 'house-prism', group: 'prism', status: 'tiles',
    name: { en: 'Pentagonal “house” prism', ru: 'Пятиугольная призма «домик»', es: 'Prisma pentagonal «casita»' },
    solid: () => prism([[0, 0], [1, 0], [1, 1], [0.5, 1 + h], [0, 1]]),
    tiling: () => ({ lattice: [[1, 0, 0], [0.5, 2 + h, 0], [0, 0, 1]], motifs: [I(), rotZ(Math.PI, [1.5, 2 + h, 0])] }),
    about: {
      en: 'The base is the equilateral pentagon “square + triangle” with angles 90°, 90°, 150°, 60°, 150°. Two adjacent angles sum to 180°, which is type 1 in the Hirschhorn–Hunt list. Rows of houses interlock roof into roof.',
      ru: 'Основание — равносторонний пятиугольник «квадрат + треугольник» с углами 90°, 90°, 150°, 60°, 150°. Два соседних угла дают 180°, это тип 1 по Хиршхорну–Ханту. Ряды «домиков» вставляются друг в друга крышами.',
      es: 'La base es el pentágono equilátero «cuadrado + triángulo» con ángulos 90°, 90°, 150°, 60°, 150°. Dos ángulos adyacentes suman 180°, el tipo 1 de la lista de Hirschhorn–Hunt. Las filas de casitas encajan tejado con tejado.',
    },
  },
  {
    id: 'gyrobifastigium', group: 'stereo', status: 'tiles',
    name: { en: 'Gyrobifastigium', ru: 'Гиробифастигиум', es: 'Girobifastigio' },
    solid: () => gyro(0),
    tiling: () => gyroTiling(0),
    about: {
      en: 'Johnson solid J26: two triangular prisms glued along a square with a quarter turn. The only Johnson solid that fills space. Layers of prisms along x and along y alternate.',
      ru: 'Многогранник Джонсона J26: две треугольные призмы, склеенные по квадрату с поворотом на 90°. Единственный многогранник Джонсона, заполняющий пространство. Слои призм вдоль x и вдоль y чередуются.',
      es: 'Sólido de Johnson J26: dos prismas triangulares pegados por un cuadrado con un giro de 90°. El único sólido de Johnson que llena el espacio. Se alternan capas de prismas a lo largo de x y de y.',
    },
  },
  {
    id: 'elongated-gyrobifastigium', group: 'stereo', status: 'tiles',
    name: { en: 'Elongated gyrobifastigium', ru: 'Удлинённый гиробифастигиум', es: 'Girobifastigio alargado' },
    solid: () => gyro(1),
    tiling: () => gyroTiling(1),
    about: {
      en: 'A cube inserted between the halves of a gyrobifastigium: the equilateral member of the “gabled rhombohedra” family (MathWorld). Same stacking, with extra layers of cubes.',
      ru: 'Между половинками гиробифастигиума вставлен куб. Равносторонний член семейства gabled rhombohedra (MathWorld). Укладка та же, только со слоями кубов.',
      es: 'Un cubo insertado entre las mitades del girobifastigio: el miembro equilátero de la familia «gabled rhombohedra» (MathWorld). El mismo apilamiento, con capas de cubos.',
    },
  },
  {
    id: 'gen-gyro', group: 'found', status: 'tiles', expect: { faces: [8, 7] },
    name: { en: 'Generalized gyrobifastigium', ru: 'Обобщённый гиробифастигиум', es: 'Girobifastigio generalizado' },
    params: [
      { key: 'phi', label: { en: 'rhombus angle φ', ru: 'угол ромба φ', es: 'ángulo del rombo φ' }, min: 40, max: 140, step: 0.5, value: 66, unit: '°' },
      { key: 'psi1', label: { en: 'upper prism tilt ψ₁', ru: 'наклон верхней призмы ψ₁', es: 'inclinación del prisma superior ψ₁' }, min: 40, max: 140, step: 0.5, value: 85, unit: '°' },
      { key: 'psi2', label: { en: 'lower prism tilt ψ₂', ru: 'наклон нижней призмы ψ₂', es: 'inclinación del prisma inferior ψ₂' }, min: 40, max: 140, step: 0.5, value: 92, unit: '°' },
    ],
    solid: (p) => genGyro(p).points,
    tiling: (p) => genGyro(p).tiling,
    about: {
      en: 'Two oblique triangular prisms (equilateral triangles, rhombic sides) glued along a common rhombus, one with rods along a, the other along b. Family types F8-203 (3 parameters) and, when a triangle becomes coplanar with a side face, F7-029. Every convex member tiles: layers of a-rods and b-rods alternate, the second tile is the point reflection x ↦ a + b + apex₁ − x, and the period between layers is apex₁ − apex₂.',
      ru: 'Две наклонные треугольные призмы (правильные треугольники, ромбические боковые грани), склеенные по общему ромбу: у одной стержни идут вдоль a, у другой вдоль b. Семейство типа F8-203 (3 параметра), а когда треугольник ложится в плоскость боковой грани — тип F7-029. Замощает каждый выпуклый член: слои a-стержней и b-стержней чередуются, вторая плитка — центральная симметрия x ↦ a + b + apex₁ − x, период между слоями apex₁ − apex₂.',
      es: 'Dos prismas triangulares oblicuos (triángulos equiláteros, caras laterales rómbicas) pegados por un rombo común, uno con barras a lo largo de a y otro a lo largo de b. Familia del tipo F8-203 (3 parámetros) y, cuando un triángulo queda coplanar con una cara lateral, del tipo F7-029. Todo miembro convexo tesela: se alternan capas de barras a y b, la segunda tesela es la simetría central x ↦ a + b + apex₁ − x y el período entre capas es apex₁ − apex₂.',
    },
  },
  {
    id: 'gen-elongated', group: 'found', status: 'tiles', expect: { faces: [8] },
    name: { en: 'Generalized elongated gyrobifastigium', ru: 'Обобщённый удлинённый гиробифастигиум', es: 'Girobifastigio alargado generalizado' },
    params: [
      { key: 'gamma', label: { en: 'angle (a, b)', ru: 'угол (a, b)', es: 'ángulo (a, b)' }, min: 50, max: 130, step: 0.5, value: 80, unit: '°' },
      { key: 'beta', label: { en: 'angle (a, c)', ru: 'угол (a, c)', es: 'ángulo (a, c)' }, min: 50, max: 130, step: 0.5, value: 75, unit: '°' },
      { key: 'alpha', label: { en: 'angle (b, c)', ru: 'угол (b, c)', es: 'ángulo (b, c)' }, min: 50, max: 130, step: 0.5, value: 70, unit: '°' },
    ],
    solid: (p) => genElongated(p).points,
    tiling: (p) => genElongated(p).tiling,
    about: {
      en: 'A rhombic parallelepiped (a, b, c) with a triangular prism on top (rods along a, triangles in the plane of b and c) and another below (rods along b, triangles in the plane of a and c). Each triangle merges with a face of the parallelepiped into a “house” pentagon: 4 pentagons + 4 rhombi, type F8-255, 3 parameters. Every member tiles with the lattice a, b, c + apex₁ − apex₂ and a point reflection.',
      ru: 'Ромбический параллелепипед (a, b, c), сверху треугольная призма (стержни вдоль a, треугольники в плоскости b и c), снизу ещё одна (стержни вдоль b, треугольники в плоскости a и c). Каждый треугольник сливается с гранью параллелепипеда в пятиугольник-«домик»: 4 пятиугольника + 4 ромба, тип F8-255, 3 параметра. Замощает каждый член: решётка a, b, c + apex₁ − apex₂ плюс центральная симметрия.',
      es: 'Un paralelepípedo rómbico (a, b, c) con un prisma triangular encima (barras a lo largo de a, triángulos en el plano de b y c) y otro debajo (barras a lo largo de b, triángulos en el plano de a y c). Cada triángulo se funde con una cara del paralelepípedo en un pentágono «casita»: 4 pentágonos + 4 rombos, tipo F8-255, 3 parámetros. Todo miembro tesela con la red a, b, c + apex₁ − apex₂ y una simetría central.',
    },
  },
  {
    id: 'rhombus-triangle', group: 'found', status: 'tiles', expect: { faces: [8] },
    name: { en: 'Rhombus ⊕ triangle', ru: 'Ромб ⊕ треугольник', es: 'Rombo ⊕ triángulo' },
    params: [
      { key: 'gamma', label: { en: 'rhombus angle', ru: 'угол ромба', es: 'ángulo del rombo' }, min: 40, max: 140, step: 0.5, value: 70, unit: '°' },
      { key: 'delta', label: { en: 'direction of the horizontal triangle side (share of the angle)', ru: 'направление горизонтальной стороны треугольника (доля угла)', es: 'dirección del lado horizontal del triángulo (fracción del ángulo)' }, min: 0.05, max: 0.95, step: 0.01, value: 0.45, unit: '' },
      { key: 'tau', label: { en: 'triangle tilt', ru: 'наклон треугольника', es: 'inclinación del triángulo' }, min: 20, max: 160, step: 0.5, value: 75, unit: '°' },
    ],
    solid: (p) => rhombusTriangle(p).points,
    tiling: (p) => rhombusTriangle(p).tiling,
    about: {
      en: 'The Minkowski sum of a unit rhombus R and an equilateral unit triangle Δ, one side of which is parallel to the plane of R. Top face R, bottom face the hexagon R ⊕ side, 4 rhombi and 2 triangles on the sides: type F8-249, 3 parameters. Every member tiles: lattice d₁, a + b + d₂, b − a, where d₁, d₂ are the lower vertices of Δ, plus the point reflection in the midpoint of an edge of R.',
      ru: 'Сумма Минковского единичного ромба R и правильного единичного треугольника Δ, одна сторона которого параллельна плоскости R. Сверху грань R, снизу шестиугольник R ⊕ сторона, по бокам 4 ромба и 2 треугольника: тип F8-249, 3 параметра. Замощает каждый член: решётка d₁, a + b + d₂, b − a, где d₁, d₂ — нижние вершины Δ, плюс центральная симметрия относительно середины ребра R.',
      es: 'La suma de Minkowski de un rombo unitario R y un triángulo equilátero unitario Δ con un lado paralelo al plano de R. Cara superior R, cara inferior el hexágono R ⊕ lado, y en los lados 4 rombos y 2 triángulos: tipo F8-249, 3 parámetros. Todo miembro tesela: red d₁, a + b + d₂, b − a, donde d₁, d₂ son los vértices inferiores de Δ, más la simetría central respecto al punto medio de una arista de R.',
    },
  },
  {
    id: 'tri-tri', group: 'found', status: 'tiles', expect: { faces: [8, 9] },
    name: { en: 'Triangle ⊕ triangle', ru: 'Треугольник ⊕ треугольник', es: 'Triángulo ⊕ triángulo' },
    params: [
      { key: 'al', label: { en: 'rotation α', ru: 'поворот α', es: 'giro α' }, min: 0, max: 180, step: 0.5, value: 40, unit: '°' },
      { key: 'be', label: { en: 'rotation β', ru: 'поворот β', es: 'giro β' }, min: 5, max: 175, step: 0.5, value: 70, unit: '°' },
      { key: 'ga', label: { en: 'rotation γ', ru: 'поворот γ', es: 'giro γ' }, min: 0, max: 180, step: 0.5, value: 25, unit: '°' },
    ],
    solid: (p) => minkowski([triangle(), rotate(triangle(), p)]),
    tiling: (p) => ({ candidates: triTriCandidates(triangle(), rotate(triangle(), p)) }),
    about: {
      en: 'The Minkowski sum of two unit equilateral triangles in general position. Depending on their relative position it has 8 faces (the generalized gyrobifastigium) or 9 faces: 4 triangles + 5 rhombi. The 9-face chamber is new to our catalog; it was found by the literature survey and confirmed by our own search. Every member tiles space with a lattice and a point reflection.',
      ru: 'Сумма Минковского двух единичных правильных треугольников общего положения. В зависимости от их взаимного расположения у неё 8 граней (обобщённый гиробифастигиум) или 9 граней: 4 треугольника и 5 ромбов. 9-гранная камера впервые попала в наш каталог; её нашёл обзор литературы, а наш собственный поиск подтвердил. Каждый член замощает пространство решёткой и центральной симметрией.',
      es: 'La suma de Minkowski de dos triángulos equiláteros unitarios en posición general. Según su posición relativa tiene 8 caras (el girobifastigio generalizado) o 9 caras: 4 triángulos y 5 rombos. La cámara de 9 caras es nueva en nuestro catálogo; la encontró la revisión bibliográfica y la confirmó nuestra propia búsqueda. Todo miembro tesela el espacio con una red y una simetría central.',
    },
  },
  {
    id: 'rhomb-tri', group: 'found', status: 'tiles', expect: { faces: [10] },
    name: { en: 'Rhombus ⊕ triangle, general position', ru: 'Ромб ⊕ треугольник, общее положение', es: 'Rombo ⊕ triángulo, posición general' },
    params: [
      { key: 'th', label: { en: 'rhombus angle', ru: 'угол ромба', es: 'ángulo del rombo' }, min: 40, max: 140, step: 0.5, value: 75, unit: '°' },
      { key: 'al', label: { en: 'rotation α', ru: 'поворот α', es: 'giro α' }, min: 0, max: 180, step: 0.5, value: 35, unit: '°' },
      { key: 'be', label: { en: 'rotation β', ru: 'поворот β', es: 'giro β' }, min: 5, max: 175, step: 0.5, value: 60, unit: '°' },
      { key: 'ga', label: { en: 'rotation γ', ru: 'поворот γ', es: 'giro γ' }, min: 0, max: 180, step: 0.5, value: 20, unit: '°' },
    ],
    solid: (p) => minkowski([rhombus(p.th), rotate(triangle(), p)]),
    tiling: (p) => ({ candidates: rhombTriCandidates(rhombus(p.th), rotate(triangle(), p)) }),
    about: {
      en: 'The Minkowski sum of a unit rhombus and a unit equilateral triangle in general position: 10 faces (2 triangles + 8 rhombi), 11 vertices, 4 parameters. It tiles space with a lattice and a point reflection. The special position with a triangle side parallel to the rhombus is the 8-face body “Rhombus ⊕ triangle” above.',
      ru: 'Сумма Минковского единичного ромба и единичного правильного треугольника общего положения: 10 граней (2 треугольника и 8 ромбов), 11 вершин, 4 параметра. Замощает пространство решёткой и центральной симметрией. Частный случай, когда сторона треугольника параллельна плоскости ромба, — 8-гранное тело «Ромб ⊕ треугольник» выше.',
      es: 'La suma de Minkowski de un rombo unitario y un triángulo equilátero unitario en posición general: 10 caras (2 triángulos y 8 rombos), 11 vértices, 4 parámetros. Tesela el espacio con una red y una simetría central. La posición especial con un lado del triángulo paralelo al rombo es el cuerpo de 8 caras «Rombo ⊕ triángulo» de arriba.',
    },
  },
  {
    id: 'tetrahedron', group: 'no', status: 'no',
    name: { en: 'Regular tetrahedron', ru: 'Правильный тетраэдр', es: 'Tetraedro regular' },
    solid: () => [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map((p) => p.map((x) => x / (2 * Math.SQRT2))),
    about: {
      en: 'The only equilateral tetrahedron. Its dihedral angle arccos(1/3) ≈ 70.53° sums to neither 360° nor 180°. This is Aristotle’s famous mistake: he believed tetrahedra fill space.',
      ru: 'Единственный равносторонний тетраэдр. Двугранный угол arccos(1/3) ≈ 70.53° не даёт в сумме ни 360°, ни 180°. Это ошибка Аристотеля, который считал, что тетраэдры заполняют пространство.',
      es: 'El único tetraedro equilátero. Su ángulo diedro arccos(1/3) ≈ 70,53° no suma ni 360° ni 180°. Es el famoso error de Aristóteles, que creía que los tetraedros llenan el espacio.',
    },
  },
  {
    id: 'square-pyramid', group: 'no', status: 'no',
    name: { en: 'Square pyramid J1', ru: 'Квадратная пирамида J1', es: 'Pirámide cuadrada J1' },
    solid: () => [[0.5, 0.5, 0], [-0.5, 0.5, 0], [-0.5, -0.5, 0], [0.5, -0.5, 0], [0, 0, Math.SQRT1_2]],
    about: {
      en: 'The only equilateral quadrilateral pyramid: the base must be a cyclic rhombus, i.e. a square. Its angles 54.74° and 109.47° = 2·54.74°, so every sum is a multiple of 54.74°, and 360°/54.74° ≈ 6.58.',
      ru: 'Единственная равносторонняя четырёхугольная пирамида: основание должно быть вписанным ромбом, то есть квадратом. Углы 54.74° и 109.47° = 2·54.74° — все суммы кратны 54.74°, а 360°/54.74° ≈ 6.58.',
      es: 'La única pirámide cuadrangular equilátera: la base debe ser un rombo inscrito, es decir, un cuadrado. Sus ángulos 54,74° y 109,47° = 2·54,74°, así que toda suma es múltiplo de 54,74°, y 360°/54,74° ≈ 6,58.',
    },
  },
  {
    id: 'octahedron', group: 'no', status: 'no',
    name: { en: 'Octahedron', ru: 'Октаэдр', es: 'Octaedro' },
    solid: () => [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]].map((p) => p.map((x) => x * Math.SQRT1_2)),
    about: {
      en: 'Dihedral angle 109.47°. Together with tetrahedra it fills space, but not on its own.',
      ru: 'Двугранный угол 109.47°. Вместе с тетраэдрами заполняет пространство, но в одиночку не может.',
      es: 'Ángulo diedro 109,47°. Junto con tetraedros llena el espacio, pero no por sí solo.',
    },
  },
  {
    id: 'triacontahedron', group: 'no', status: 'no',
    name: { en: 'Rhombic triacontahedron', ru: 'Ромботриаконтаэдр', es: 'Triacontaedro rómbico' },
    zonohedron: () => [[0, 1, phi], [0, -1, phi], [1, phi, 0], [-1, phi, 0], [phi, 0, 1], [phi, 0, -1]].map(unit),
    about: {
      en: '30 golden rhombi. All dihedral angles equal 144°, which divides neither 360° nor 180°.',
      ru: '30 золотых ромбов. Все двугранные углы равны 144°, а 144° не делит ни 360°, ни 180°.',
      es: '30 rombos áureos. Todos los ángulos diedros miden 144°, que no divide ni 360° ni 180°.',
    },
  },
  {
    id: 'engel38', group: 'records', status: 'tiles', ...stereo('engel38'),
    name: { en: 'Engel’s 38-hedron (record)', ru: '38-гранник Энгеля (рекорд)', es: 'El 38-edro de Engel (récord)' },
    about: {
      en: 'The convex space-filler with the most faces known: 38 faces, 70 vertices, 106 edges (P. Engel, 1981). It is the Voronoi cell of a point in the cubic space group I4₁32; its 48 symmetric copies per unit cell fill space. We reproduced it exactly (rational certificate) and proved that it is rigid: no nearby non-Voronoi tile with the same structure exists. Whether 38 can be beaten is open.',
      ru: 'Выпуклый многогранник, заполняющий пространство, с наибольшим известным числом граней: 38 граней, 70 вершин, 106 рёбер (П. Энгель, 1981). Это ячейка Вороного точки в кубической группе I4₁32; её 48 симметричных копий на ячейку решётки заполняют пространство. Мы точно воспроизвели его (рациональный сертификат) и доказали жёсткость: близкой не-вороновской плитки с той же структурой нет. Можно ли превзойти 38 — открытый вопрос.',
      es: 'El poliedro convexo que llena el espacio con más caras conocido: 38 caras, 70 vértices, 106 aristas (P. Engel, 1981). Es la celda de Voronói de un punto en el grupo cúbico I4₁32; sus 48 copias simétricas por celda llenan el espacio. Lo reprodujimos exactamente (certificado racional) y probamos que es rígido: no existe una tesela no Voronói cercana con la misma estructura. Si 38 puede superarse sigue abierto.',
    },
  },
  {
    id: 'schmitt35', group: 'records', status: 'tiles', ...stereo('schmitt35'),
    name: { en: '35-hedron in I4₁22', ru: '35-гранник в I4₁22', es: '35-edro en I4₁22' },
    about: {
      en: 'The best tetragonal example (M. Schmitt, 2016): a Voronoi cell with 35 faces in the space group I4₁22 at c/a = 1.454. Reproduced exactly by our code.',
      ru: 'Лучший тетрагональный пример (М. Шмитт, 2016): ячейка Вороного с 35 гранями в группе I4₁22 при c/a = 1,454. Точно воспроизведён нашим кодом.',
      es: 'El mejor ejemplo tetragonal (M. Schmitt, 2016): una celda de Voronói con 35 caras en el grupo I4₁22 con c/a = 1,454. Reproducido exactamente por nuestro código.',
    },
  },
  {
    id: 'ours33', group: 'records', status: 'tiles', ...stereo('ours33'),
    name: { en: '33-hedron in I-42d (new)', ru: '33-гранник в I-42d (новый)', es: '33-edro en I-42d (nuevo)' },
    about: {
      en: 'Found by our search: a Voronoi cell with 33 faces in the group I-42d; the previous best for this group was 31 (Schmitt 2016). Certified in exact rational arithmetic.',
      ru: 'Найден нашим поиском: ячейка Вороного с 33 гранями в группе I-42d; прежний лучший результат для этой группы — 31 (Шмитт, 2016). Подтверждён точной рациональной арифметикой.',
      es: 'Encontrado por nuestra búsqueda: una celda de Voronói con 33 caras en el grupo I-42d; el mejor anterior para este grupo era 31 (Schmitt, 2016). Certificado con aritmética racional exacta.',
    },
  },
  {
    id: 'ours29', group: 'records', status: 'tiles', ...stereo('ours29'),
    name: { en: '29-hedron in I4₁/a (new)', ru: '29-гранник в I4₁/a (новый)', es: '29-edro en I4₁/a (nuevo)' },
    about: {
      en: 'Found by our search: a Voronoi cell with 29 faces in the group I4₁/a; the previous best for this group was 23 (Schmitt 2016). Certified in exact rational arithmetic.',
      ru: 'Найден нашим поиском: ячейка Вороного с 29 гранями в группе I4₁/a; прежний лучший результат для этой группы — 23 (Шмитт, 2016). Подтверждён точной рациональной арифметикой.',
      es: 'Encontrado por nuestra búsqueda: una celda de Voronói con 29 caras en el grupo I4₁/a; el mejor anterior para este grupo era 23 (Schmitt, 2016). Certificado con aritmética racional exacta.',
    },
  },
];

// ---------- семейства, найденные перебором ----------
const V3 = {
  add: (...vs) => vs.reduce((s, x) => [s[0] + x[0], s[1] + x[1], s[2] + x[2]], [0, 0, 0]),
  mul: (x, k) => x.map((c) => c * k),
  sub: (x, y) => [x[0] - y[0], x[1] - y[1], x[2] - y[2]],
  dot: (x, y) => x[0] * y[0] + x[1] * y[1] + x[2] * y[2],
  cross: (x, y) => [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]],
  unit: (x) => { const l = Math.hypot(...x); return x.map((c) => c / l); },
};
const reflect = (w) => new THREE.Matrix4().makeScale(-1, -1, -1).setPosition(...w);

function genGyro({ phi, psi1, psi2 }) {
  const a = [1, 0, 0], b = [Math.cos(rad(phi)), Math.sin(rad(phi)), 0], z = [0, 0, 1];
  const mb = V3.cross(z, b), ma = V3.cross(z, a);
  const u1 = V3.add(V3.mul(mb, Math.cos(rad(psi1))), V3.mul(z, Math.sin(rad(psi1))));
  const u2 = V3.sub(V3.mul(ma, Math.cos(rad(psi2))), V3.mul(z, Math.sin(rad(psi2))));
  const apex1 = V3.add(V3.mul(b, 0.5), V3.mul(u1, h));
  const apex2 = V3.add(V3.mul(a, 0.5), V3.mul(u2, h));
  const O = [0, 0, 0];
  const T1 = [O, b, apex1], T2 = [O, a, apex2];
  return {
    points: [...T1, ...T1.map((p) => V3.add(p, a)), ...T2, ...T2.map((p) => V3.add(p, b))],
    tiling: { lattice: [a, b, V3.sub(apex1, apex2)], motifs: [I(), reflect(V3.add(a, b, apex1))] },
  };
}

function genElongated({ gamma, beta, alpha }) {
  const a = [1, 0, 0], b = [Math.cos(rad(gamma)), Math.sin(rad(gamma)), 0];
  const cx = Math.cos(rad(beta)), cy = (Math.cos(rad(alpha)) - cx * b[0]) / b[1];
  const c = [cx, cy, Math.sqrt(Math.max(1e-6, 1 - cx * cx - cy * cy))];
  const u1 = V3.unit(V3.sub(c, V3.mul(b, V3.dot(c, b))));
  const u2 = V3.unit(V3.sub(c, V3.mul(a, V3.dot(c, a))));
  const apex1 = V3.add(c, V3.mul(b, 0.5), V3.mul(u1, h));
  const apex2 = V3.sub(V3.mul(a, 0.5), V3.mul(u2, h));
  const O = [0, 0, 0];
  const box = [O, a, b, V3.add(a, b), c, V3.add(a, c), V3.add(b, c), V3.add(a, b, c)];
  return {
    points: [...box, apex1, V3.add(apex1, a), apex2, V3.add(apex2, b)],
    tiling: { lattice: [a, b, V3.add(c, V3.sub(apex1, apex2))], motifs: [I(), reflect(V3.add(a, b, c, apex1))] },
  };
}

// ---------- суммы Минковского ----------
function triangle() {
  const R = 1 / Math.sqrt(3);
  return [0, 1, 2].map((k) => [R * Math.cos((2 * Math.PI * k) / 3), R * Math.sin((2 * Math.PI * k) / 3), 0]);
}

function rhombus(th) {
  const c = Math.cos(rad(th)), s = Math.sin(rad(th));
  return [[0, 0, 0], [1, 0, 0], [1 + c, s, 0], [c, s, 0]].map((v) => [v[0] - (1 + c) / 2, v[1] - s / 2, 0]);
}

function rotate(P, { al, be, ga }) {
  const [a, b, g] = [rad(al), rad(be), rad(ga)];
  const Rx = (v) => [v[0], Math.cos(a) * v[1] - Math.sin(a) * v[2], Math.sin(a) * v[1] + Math.cos(a) * v[2]];
  const Ry = (v) => [Math.cos(b) * v[0] + Math.sin(b) * v[2], v[1], -Math.sin(b) * v[0] + Math.cos(b) * v[2]];
  const Rz = (v) => [Math.cos(g) * v[0] - Math.sin(g) * v[1], Math.sin(g) * v[0] + Math.cos(g) * v[1], v[2]];
  return P.map((v) => Rz(Ry(Rx(v))));
}

function minkowski(parts) {
  return parts.reduce((acc, P) => acc.flatMap((x) => P.map((y) => V3.add(x, y))));
}

function triLabels(T) {
  const out = [];
  for (let i = 0; i < 3; i++) {
    const o = T[i];
    const rest = [0, 1, 2].filter((j) => j !== i).map((j) => V3.sub(T[j], o));
    out.push([o, rest[0], rest[1]], [o, rest[1], rest[0]]);
  }
  return out;
}

// Разбиения суммы двух треугольников: камера 8 граней — формула A, камера 9 граней — формула B.
function triTriCandidates(T1, T2) {
  const out = [];
  for (const [o1, u1, u2] of triLabels(T1)) {
    for (const [o2, v1, v2] of triLabels(T2)) {
      const c = V3.mul(V3.add(o1, o2), 2);
      out.push({ lattice: [v1, u1, V3.sub(u2, v2)], motifs: [I(), reflect(V3.add(c, v1, u1, u2))] });
      out.push({ lattice: [V3.add(u1, v1), V3.add(u1, v2), V3.sub(u2, u1)], motifs: [I(), reflect(c)] });
    }
  }
  return out;
}

function rhombTriCandidates(R, T) {
  const out = [];
  for (let i = 0; i < 4; i++) {
    const o = R[i];
    const nb = [V3.sub(R[(i + 1) % 4], o), V3.sub(R[(i + 3) % 4], o)];
    for (const [r1, r2] of [nb, [nb[1], nb[0]]]) {
      for (const [ot, t1, t2] of triLabels(T)) {
        const r12 = V3.add(r1, r2);
        out.push({
          lattice: [V3.add(r12, t1), V3.add(r12, t2), V3.add(V3.sub(r1, r2), V3.sub(t1, t2))],
          motifs: [I(), reflect(V3.add(V3.mul(V3.add(o, ot), 2), r1, t1))],
        });
      }
    }
  }
  return out;
}

function rhombusTriangle({ gamma, delta, tau }) {
  const a = [1, 0, 0], b = [Math.cos(rad(gamma)), Math.sin(rad(gamma)), 0], z = [0, 0, 1];
  const dir = rad(gamma * delta);
  const c = [Math.cos(dir), Math.sin(dir), 0];
  const m = V3.cross(z, c);
  const w = V3.sub(V3.mul(m, Math.cos(rad(tau))), V3.mul(z, Math.sin(rad(tau))));
  const d1 = V3.add(V3.mul(c, -0.5), V3.mul(w, h));
  const d2 = V3.add(d1, c);
  const O = [0, 0, 0];
  const R = [O, a, V3.add(a, b), b];
  const points = R.flatMap((r) => [r, V3.add(r, d1), V3.add(r, d2)]);
  return { points, tiling: { lattice: [d1, V3.add(a, b, d2), V3.sub(b, a)], motifs: [I(), reflect(a)] } };
}

function gyro(e) {
  const z = e / 2;
  const sq = (zz) => [[0.5, 0.5, zz], [-0.5, 0.5, zz], [-0.5, -0.5, zz], [0.5, -0.5, zz]];
  return [
    ...sq(z), ...(e ? sq(-z) : []),
    [0.5, 0, z + h], [-0.5, 0, z + h],
    [0, 0.5, -z - h], [0, -0.5, -z - h],
  ];
}

function gyroTiling(e) {
  return {
    lattice: [[1, 0, 0], [0, 1, 0], [0.5, 0, 2 * (e + h)]],
    motifs: [I(), rotZ(Math.PI / 2, [0, 0.5, e + h])],
  };
}
