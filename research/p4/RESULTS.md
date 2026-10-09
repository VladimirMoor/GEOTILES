# Problem 4: a convex space-filler with more than 38 faces

## Literature check (2026-10-09)

**Record: 38 faces.** P. Engel, *Z. Kristallogr.* 154 (1981) 199–215, found a Dirichlet–Voronoi (DV) stereohedron, i.e. a plesiohedron, for the cubic group I4₁32 (IT 214). Its f-vector is (70, 106, 38). There are four combinatorial types with these numbers. 38 is the largest number of faces known for any convex monohedral tile in ℝ³.

**What is actually proved.** Wikipedia and MathWorld say that Schmitt (2016) proved 38 to be the maximum for plesiohedra. **That is not what the thesis says.**

M. Schmitt, *On Space Groups and Dirichlet–Voronoi Stereohedra*, PhD thesis, FU Berlin 2016, §2.3:
- It is a **numerical grid search** over the fundamental domains of the normalisers for 145 space groups: the tetragonal, trigonal, hexagonal and cubic ones. Triclinic, monoclinic and orthorhombic groups are omitted; they cannot produce more than 38 faces.
- It finds 3315 combinatorial types of DV-stereohedra. Every facet number from 4 to 38 occurs, and the maximum is attained only by Engel's examples.
- **It contains no proof that 38 is the maximum.**
- Schmitt himself calls problem (G) of Grünbaum–Shephard (1980) open: the least upper bound for the number of facets of a convex monohedral prototile in ℝⁿ. For n ≥ 3 it is not even known that the number is bounded.

**Rigorous upper bounds for DV-stereohedra**, proved group by group (Bochiş–Santos 2001, 2006; Sabariego–Santos 2008, 2011):

| Groups | Bound |
|---|---|
| Groups with reflections | ≤ 18 |
| Non-cubic groups | ≤ 80 |
| Full cubic groups | ≤ 25 |
| Quarter cubic groups | ≤ 92 (arXiv:0708.2114) |

So for plesiohedra the true maximum lies somewhere between 38 and 92.

**General bounds.**
- Every stereohedron (the tile of an isohedral tiling): Delone's bound 390, improved by Tarasov to 378.
- Non-isohedral monohedral tilings (anisohedral, k-isohedral): **no bound is known.**

**Not a record.** The 2026 "Josehedron" (arXiv:2604.07160) is a 12-face plesiohedron.

## Where room is left

1. **Plesiohedra, i.e. DV cells of one orbit.**
   - Schmitt's grid search covered this class thoroughly, numerically.
   - Possible gaps: very small parameter chambers that a grid misses, and special positions.
   - Low expected yield.
2. **Stereohedra that are not Dirichlet cells.**
   - These are convex tiles of isohedral tilings whose shape is not the Voronoi cell of an orbit. The gyrobifastigium is an example.
   - Bound 378; we know of **no systematic search**.
3. **Non-isohedral (k-isohedral) convex monohedral tilings.**
   - Example: Voronoi cells of several orbits that happen to be congruent.
   - **No bound and no systematic search known.**

## Search 1: DV-stereohedra, all groups 75–230 (2026-10-10)

**Tools.**
- `dvfast.py`: DV cells with a fixed, provably sufficient neighbour radius, namely twice the covering radius of the lattice.
- `exact.py`: an exact rational certificate that a cell has at least k facets. For each facet it exhibits three non-collinear exactly feasible points lying on that facet's plane. It uses fractional coordinates with the metric tensor, so hexagonal cells are handled too.
- `batch.py`: memory-safe batch runs.
- `verify_batch.py`: exact re-check of the best points for every group.

**Validation.** Schmitt's reference points are reproduced exactly:
- IT 214: (70, 106, 38);
- IT 98 at c/a = 1.454: 35.

**Run.** For each of the 156 groups:
- 8 000 random points, with c/a ∈ [0.3, 4] for non-cubic groups;
- then 25 local hill-climbs × 700 steps that steer towards chamber walls;
- then exact certification of the best points.

**Warning about numerical artefacts.** Float maxima of 58, 44 and 38 for IT 141, 166 and 134 were artefacts. The points sat about 10⁻¹² from a mirror plane, and nearly coincident planes were counted as separate facets. The exact values are 29, 22 and 19, exactly Schmitt's values. **Only certified numbers count.**

**Results** (full table: `verified_batch1.txt`).
- **No cell with more than 38 facets.**
- Certified maxima **equal** to Schmitt's: IT 214 (38), 178 (34), 141 (29), 206 (28), 166 (22), 134 (19), 121 (17), and others.
- Certified maxima **above** Schmitt's, i.e. new lower bounds for these groups. Schmitt's values were re-checked against the thesis text.

  | IT | group | certified | Schmitt | c/a of our point |
  |---|---|---|---|---|
  | 86 | P4₂/n | 21 | 20 | 0.588 |
  | 88 | I4₁/a | **29** | 23 | 0.606 |
  | 122 | I-42d | 33 | 31 | 1.487 |
  | 167 | R-3c | 27 | 26 | 2.229 |
  | 201 | Pn-3 | 17 | 16 | cubic |
  | 228 | Fd-3c | 17 | 16 | cubic |
  | 138 | P4₂/ncm | 17 | 14 | 0.425, outside Schmitt's range [0.5, 3.5] |
  | 184 | P6cc | 13 | 11 | 0.300, at the edge of our range |
  | 192 | P6/mcc | 13 | 11 | 0.328, outside Schmitt's range |

- In 56 groups our short search stayed below Schmitt. The search per group is small, and Schmitt's grid is about 10⁹ points.

**Next.** Schmitt sampled c/a only in [1/2, 7/2]. The proven bound for non-cubic groups is 80. So we search the extreme metrics, c/a ∈ [0.05, 0.5] ∪ [3.5, 20], for groups without reflections.
