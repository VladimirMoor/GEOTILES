# Задача 1. Равносторонние выпуклые многогранники, замощающие пространство: обзор литературы

Дата поиска: 2026-10-07.

## Вывод

Классификации равносторонних (все рёбра равны) выпуклых многогранников, замощающих R^3, в литературе не найдено.
Двумерный аналог закрыт: Hirschhorn–Hunt (1985), упрощённое доказательство у Klaassen (arXiv:2506.18473, 2025, на основе Rao 2017).
Klaassen о 3D не пишет. В обзоре 2026 года (arXiv:2604.18545) прямо сказано, что общей классификации нет даже для тетраэдров.

## Что известно (частичные результаты)

| Класс | Результат | Источник |
|---|---|---|
| Правильногранные | Ровно 5 многогранников заполняют пространство: куб, треуг. призма, шестиуг. призма, усечённый октаэдр, гиробифастигиум (J26) | Johnson 1966, классификация правильногранных |
| Параллелоэдры (Фёдоров, 5 типов) | У каждого типа есть равносторонняя реализация (равносторонний зоноэдр: все образующие одной длины). Ромбоэдры, ромбододекаэдры Кеплера и Билинского, удлинённый додекаэдр, усечённый октаэдр | Fedorov 1885; Bilinski 1960; Grünbaum 2010 |
| Призмы | Призма с ромбическими боковыми гранями над равносторонним многоугольником, замощающим плоскость → заполняет пространство. Основания классифицированы: треугольник, ромб, пятиугольники Hirschhorn–Hunt, шестиугольники Klaassen | Hirschhorn–Hunt 1985; Klaassen 2025 |
| Gabled rhombohedron | Семейство удлинённых гиробифастигиумов, есть равносторонний член | MathWorld |
| Каталоги Голдберга по числу граней | 5–10 граней (не равносторонние; это списки известных примеров, полнота не доказана) | Goldberg, JCT 1972/74; Geom. Dedicata 1977, 1978, 1981, 1982 |
| Изоэдральные простые разбиения | Комбинаторно перечислены все для плиток с ≤16 гранями; плиток с <14 гранями нет | Delgado-Friedrichs, O'Keeffe, Acta Cryst A 2005 |

## Наблюдения для постановки

- Пространство реализаций 3-многогранника с e рёбрами по модулю подобия имеет размерность e−1, а равносторонность даёт e−1 уравнений. Поэтому для большинства комбинаторных типов равносторонних реализаций конечное число, и проверка замощения сводится к конечному числу конкретных тел. Исключения — гибкие типы: ромбические грани, призмы.
- Сильный необходимый фильтр: вокруг каждого ребра тайлинга двугранные углы плиток в сумме дают 2π (или π, если ребро лежит внутри грани соседа).
- Малые случаи сразу:
  - n=4: только правильный тетраэдр, не замощает.
  - n=5: правильная квадратная пирамида J1 (не замощает, углы 54.74°/109.47° не складываются в 180°/360°) и все равносторонние треугольные призмы (замощают).
  - n=6: в частности все ромбоэдрические параллелепипеды (замощают); остальные комбинаторные типы шестигранников нужно разбирать.

## Ссылки

- Klaassen, Which equilateral convex polygons tile the plane? https://arxiv.org/abs/2506.18473
- Softening locally polyhedral tilings (2026). https://arxiv.org/pdf/2604.18545
- Space-filling polyhedron. https://en.wikipedia.org/wiki/Space-filling_polyhedron
- Gyrobifastigium. https://en.wikipedia.org/wiki/Gyrobifastigium
- Bilinski dodecahedron. https://en.wikipedia.org/wiki/Bilinski_dodecahedron
- Gabled rhombohedron. https://mathworld.wolfram.com/GabledRhombohedron.html
- Goldberg, On the space-filling hexahedra. https://link.springer.com/article/10.1007/BF00181585
- Goldberg, On the space-filling heptahedra. https://link.springer.com/article/10.1007/BF00181630
- Goldberg, On the space-filling octahedra. https://link.springer.com/article/10.1007/BF01447431
- Goldberg, On the space-filling enneahedra. https://link.springer.com/article/10.1007/BF00147314
- Delgado-Friedrichs, O'Keeffe, Isohedral simple tilings: binodal and by tiles with ≤16 faces. https://researchportalplus.anu.edu.au/en/publications/isohedral-simple-tilings-binodal-and-by-tiles-with-16-faces/
- Schein, Gayed, Convex equilateral polyhedra with polyhedral symmetry (PNAS 2014). https://www.pnas.org/doi/10.1073/pnas.1310939111
