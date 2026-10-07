// Каталог равносторонних выпуклых многогранников (длина ребра 1). Тексты на трёх языках.
// zonohedron(p) → образующие, или solid(p) → вершины; tiling(p) → 'parallelohedron' | { lattice, motifs } | нет.
import * as THREE from 'three';

const h = Math.sqrt(3) / 2;
const phi = (1 + Math.sqrt(5)) / 2;
const rad = (d) => (d * Math.PI) / 180;
const unit = (a) => { const l = Math.hypot(...a); return a.map((x) => x / l); };
const prism = (base, H = 1) => [...base.map(([x, y]) => [x, y, 0]), ...base.map(([x, y]) => [x, y, H])];
const rotZ = (a, t = [0, 0, 0]) => new THREE.Matrix4().makeRotationZ(a).setPosition(...t);
const I = () => new THREE.Matrix4();

export const GROUPS = [
  { id: 'par', title: { en: 'Parallelohedra', ru: 'Параллелоэдры', es: 'Paraleloedros' },
    note: { en: 'Tile space by translations alone', ru: 'Замощают пространство параллельными переносами', es: 'Teselan el espacio solo por traslaciones' } },
  { id: 'prism', title: { en: 'Prisms', ru: 'Призмы', es: 'Prismas' },
    note: { en: 'The base tiles the plane; layers stack up', ru: 'Основание замощает плоскость, слои укладываются стопкой', es: 'La base tesela el plano y las capas se apilan' } },
  { id: 'stereo', title: { en: 'Other space-fillers', ru: 'Другие замощающие', es: 'Otros que llenan el espacio' },
    note: { en: 'Rotations are needed', ru: 'Нужны повороты', es: 'Se necesitan rotaciones' } },
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
];

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
