# Survey: monohedral space-fillers in nature and in the literature

2026-10-08. Question: do crystals, chemistry, foams or biology already provide convex polyhedra that tile ℝ³ with congruent copies, and are any of them equilateral? Special focus: an equilateral convex space-filler with ≤ 8 faces missing from our enumeration (`RESULTS.md`), and candidates with 9–14 faces.

Method: web/literature search, plus small computations. The scripts live in the session scratchpad and are not committed. They compute Voronoi cells of single-orbit crystal structures, face counts of Minkowski sums, and exact overlap-volume checks of lattice + point-reflection tilings.

## 0. Bottom line

1. **Nothing is missing at ≤ 8 faces.** Every ≤ 8-face convex monohedral tile from natural or classical sources is either non-equilateral, or is equilateral and falls into one of our families (§7).
2. **Three of our families are Minkowski sums of unit polygons** (§8). This was established numerically here; it has not been proved.
   - F8-203, the generalized gyrobifastigium, is triangle ⊕ triangle.
   - F7-029 is the same sum in special position.
   - F8-249 is rhombus ⊕ triangle in special position.
   - Pentagonal prisms whose base has two *adjacent* angles summing to 180° are rhombus ⊕ triangle with a rhombus side parallel to the triangle plane.
3. **The Schmitt–Conway–Danzer (SCD) aperiodic biprism has a 1-parameter equilateral subfamily inside F8-203.** That subfamily needs λ = ½ and c = (√3/2)·a. Our F8-203 tilings use a point reflection, which is orientation-reversing, so they are periodic. This is consistent with SCD aperiodicity, which holds only when mirror images are excluded (§5).
4. **New equilateral space-fillers with 9 and 10 faces** (§8), found by computation here:
   - Generic triangle ⊕ triangle: 9 faces (4 T + 5 R), 3 parameters.
   - Generic rhombus ⊕ triangle: 10 faces (2 T + 8 R), 4 parameters.

   Both tile by a lattice plus a point reflection. An exact overlap-volume check confirmed this for every tested member: 14 members of the 9-face family and 12 of the 10-face family.
5. **Hexagonal prism caveat.** The short list in the task says "centrally symmetric hexagon base", but equilateral hexagons satisfying Reinhardt type 2 (A + B + D = 360°) are generally *not* centrally symmetric. They tile the plane, so the right prisms over them tile space (§7). `RESULTS.md` already covers this case through Klaassen's criterion. This note only flags the summary wording.
6. **Negative results.** The hcp Voronoi cell (trapezo-rhombic dodecahedron) and the honeybee cell have **no** equilateral convex realization in their combinatorial type (proofs in §10).
7. **Monohedral natural examples are rare.** Foams, clathrates, Frank–Kasper phases and most zeolites need two or more tile types. The genuinely monohedral natural examples are:
   - the Voronoi or Wigner–Seitz cells of single-orbit structures (§1, §2);
   - the Kelvin and Williams foams;
   - the sodalite (SOD) cage tiling.

## 1. Wigner–Seitz cells of the 14 Bravais lattices

Every Wigner–Seitz (WS) cell is one of the 5 combinatorial parallelohedra (Fedorov), and all of them tile by translations.

| Lattice(s) | WS cell | F | Equilateral? |
|---|---|---|---|
| sc, primitive tetragonal and orthorhombic, triclinic (some cases) | parallelepiped | 6 | sc: cube, yes. Rhombic parallelepipeds are F6-007 |
| hexagonal; C-centred orthorhombic; primitive monoclinic | hexagonal prism | 8 | regular prism with h = a: yes (F8-192) |
| fcc; rhombohedral (one range of α); body-centred orthorhombic with c² = a² + b² | rhombic dodecahedron (RD) | 12 | fcc: yes. Any zonotope of 4 unit vectors in general position is an equilateral RD-type parallelohedron (5 parameters), Bilinski dodecahedron included |
| body-centred tetragonal with c/a > √2; body-centred orthorhombic with c² > a² + b² | elongated dodecahedron | 12 | equilateral version exists (unit elongation of an equilateral RD) |
| bcc; body-centred tetragonal with c/a < √2; face-centred orthorhombic; rhombohedral (other range) | truncated octahedron (TO) | 14 | bcc: yes (Archimedean). bct with c/a = 1.076 (In) has 2 edge lengths |

Source: [Wikipedia, Wigner–Seitz cell](https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell) (table of the 24 Voronoi sorts in 5 classes).

## 2. Voronoi cells of single-orbit crystal structures (computed)

When all atoms form one orbit of the space group, the Voronoi cells are congruent and give a monohedral, isohedral tiling (a Dirichlet stereohedron, or plesiohedron). Cells were computed with standard structure parameters. "Edges/NN" is the range of edge lengths relative to the nearest-neighbour distance.

| Structure (examples) | F | V | E | Face types | Edges/NN | Distinct lengths | Equilateral? |
|---|---|---|---|---|---|---|---|
| sc (α-Po) | 6 | 8 | 12 | 6 squares | 1.000 | 1 | yes, cube |
| **pyrochlore lattice** (Cu in MgCu₂, 16c of Fd-3m) | 6 | 8 | 12 | 6 rhombi | 1.225 | 1 | **yes**: a rhombohedron (D3d) in 4 orientations, so the tiling is not a lattice tiling. Already in F6-007 |
| fcc (Cu) | 12 | 14 | 24 | 12 rhombi | 0.612 | 1 | yes (RD) |
| hcp, ideal c/a (Mg) | 12 | 14 | 24 | 6 rhombi + 6 trapezoids | 0.408–0.816 | 3 (ratio 2:3:4) | no; impossible in this combinatorial type (§10) |
| α-U (A20, Cmcm 4c) | 12 | 18 | 28 | 8 quads + 4 hexagons | 0.42–0.72 | 5 | no. Same f-vector and vertex degrees (16×3, 2×4) as the elongated dodecahedron, but **not centrally symmetric**: a non-lattice ED-like tiler. Candidate (§9) |
| lonsdaleite, ideal (ice Ih O-sublattice) | 11 | 15 | 24 | 7 triangles + 3 hexagons + 1 nonagon (coplanar merge at the ideal geometry) | 0.47–0.82 | 3 | no |
| bcc (W, α-Fe) | 14 | 24 | 36 | 6 squares + 8 hexagons | 0.408 | 1 | yes (TO) |
| bct, c/a = 1.076 (In, A6) | 14 | 24 | 36 | 6 + 8 | 0.335–0.445 | 2 | no; affine image of TO |
| diamond (C, Si, ice Ic) | 16 | 16 | 30 | 4 hexagons + 12 isosceles triangles | 0.50–0.82 | 2 | no; equal edges force the regular tetrahedron (§10) |
| Se/Te (A8) | 16 | 26 | 40 | 3⁴4⁶7⁴8² | 0.11–1.35 | 8 | no |
| α-Ga (A11) | 16 | 28 | 42 | 4⁸6⁷10¹ | 0.16–1.08 | 8 | no |
| srs / Laves graph (I4₁32 8a) | 17 | 30 | 45 | 4⁶5⁶6²8³ | 0.39–1.22 | 3 | no; matches the known 17-face plesiohedron |
| β-Sn (A5) | 18 | 32 | 48 | 3⁸4²6⁴10⁴ | 0.17–0.75 | 3 | no |
| α-As (A7) | 19 | 31 | 48 | 3¹²6¹8³10³ | 0.02–1.01 | 5 | no |

Literature checks: the diamond cell (triakis truncated tetrahedron, 16 faces), the hcp cell (trapezo-rhombic dodecahedron, 12 faces) and the Laves-graph cell (17 faces) all agree with [Plesiohedron](https://en.wikipedia.org/wiki/Plesiohedron), [Trapezo-rhombic dodecahedron](https://en.wikipedia.org/wiki/Trapezo-rhombic_dodecahedron), [Triakis truncated tetrahedron](https://en.wikipedia.org/wiki/Triakis_truncated_tetrahedron) and [Schoen, double Laves Voronoi cell](https://www.schoengeometry.com/e70-tpms-media/double_Laves_voronoy_cell.pdf).

Delaunay side: the bcc Delaunay tiling is monohedral, made of tetragonal disphenoids (Sommerville No. 1) with edges a, a, (√3/2)a ×4. It is not equilateral. The only equilateral tetrahedron is regular, and it does not tile.

## 3. Nets and tilings (RCSR, EPINET, ToposPro)

- **Isohedral simple tilings** are face-to-face tilings in which every edge lies in 3 tiles and every vertex in 4. [Delgado-Friedrichs & O'Keeffe, Acta Cryst. A61 (2005) 358](https://doi.org/10.1107/S0108767305009578) state that "there are no such tilings by polyhedra with less than 14 faces". At 14, 15 and 16 faces there are 10, 65 and 434 simple space-filling polyhedra, giving 23, 136 and 710 tilings. The TO and, presumably, the Williams β-tetrakaidecahedron are among the ten 14-face tiles.
  - Consequence for us: any isohedral tiling by a tile with ≤ 13 faces must be non-simple. It has vertices shared by more than 4 tiles, or it is not face-to-face. This holds for all our tilers.
  - The ten 14-face simple isohedral tiles are a ready-made list to test for equilateral realizations.
- **Tile-transitive tilings by (topological) tetrahedra**: [Delgado-Friedrichs & Huson, DCG 2000](https://link.springer.com/article/10.1007/s004540010035) found 149 types, 9 of them by combinatorial regular tetrahedra. These are combinatorial, not convex-geometric. A convex equilateral tetrahedron is regular, so this gives nothing for us.
- **Transitivity [1 1 1 1]** (vertex-, edge-, face- and tile-transitive) natural tilings: 7 examples ([Delgado-Friedrichs, O'Keeffe, Proserpio, Treacy, Acta Cryst. A79 (2023) 192](https://doi.org/10.1107/S2053273323000414)). The convex ones among these are cube-like. Most natural tiles (dia, srs) have non-planar rings and are not convex polyhedra.
- [RCSR](https://rcsr.anu.edu.au) lists tilings and their transitivities. Natural tiles are combinatorial objects, and their convex straight-face realizations usually coincide with the parallelohedra or the stereohedra above.

## 4. Zeolites, clathrates, Frank–Kasper phases, foams

| System | Tiles | Monohedral? | Equilateral? |
|---|---|---|---|
| Sodalite **SOD** (natural tile *sod*, [4⁶6⁸]) | truncated octahedron | **yes** | yes |
| LTA (*lta*, *sod*, *d4r*), FAU (*fau*, *sod*, *d6r*), CHA (*cha*, *d6r* = hexagonal prism), RHO (*lta*, *d8r*) | 2–3 tile types | no | — |
| Other zeolites ([Anurova et al., J. Phys. Chem. C 2010](https://doi.org/10.1021/jp1030027), 194 frameworks) | natural building units | we found no other tile-transitive case (not verified per framework) | — |
| Clathrate I = Weaire–Phelan = A15 (β-W) | 5¹² + 5¹²6² (Z12 + Z14) | no | — |
| Clathrate II (C15) | 5¹² + 5¹²6⁴ | no | — |
| Clathrate H | 5¹², 4³5⁶6³, 5¹²6⁸ | no | — |
| Frank–Kasper phases (A15, C14, C15, σ, μ, …) | Z12, Z14, Z15, Z16 | never monohedral | — |
| **Kelvin foam** (bcc) | flat-faced TO | **yes** | yes |
| **Williams foam** (tetragonal, 2 cells per unit cell) | β-tetrakaidecahedron: 2 squares + 8 pentagons + 4 hexagons, 24 vertices | **yes** | unknown. Candidate (§9) |
| fcc and hcp "foams" (not Plateau-stable) | RD; trapezo-rhombic dodecahedron | yes | RD yes; TRD impossible |

Sources: [Wikipedia, A15 phases](https://en.wikipedia.org/wiki/A15_phases); [Williams foam elastic constants, DTU](https://orbit.dtu.dk/en/publications/space-filling-polyhedra-as-mechanical-models-for-solidified-dry-f-2/); [Robert Williams (geometer)](https://en.wikipedia.org/wiki/Robert_Williams_(geometer)).

## 5. Stereohedra, classical catalogues, special tiles

| Item | F | Notes | Equilateral? |
|---|---|---|---|
| **SCD biprism** (Schmitt 1988; Conway–Danzer) | 8 | conv(0, a, b, a+b, c, a+c, d, b+d) with c = λb + (0,0,c), d = λa − (0,0,c), \|a\| = \|b\|. Two triangular prisms on a rhombus with angle π/2 − φ ([Baake & Frettlöh](https://arxiv.org/abs/math-ph/0411052)) | equilateral iff λ = ½ and c = (√3/2)\|a\|, a 1-parameter family (rhombus angle) **inside F8-203**. With proper motions only it tiles only aperiodically when φ ∉ πℚ ([MathWorld](https://mathworld.wolfram.com/Space-FillingPolyhedron.html): "mirror images are excluded"). Our F8-203 tiling uses a point reflection, so the two results are consistent |
| Gabled rhombohedron, Goldberg type 8-VI (= elongated gyrobifastigium) | 8 | 4 pentagons + 4 quadrilaterals, all vertices of degree 3 | yes: this is F8-255 ([Wikipedia](https://en.wikipedia.org/wiki/Elongated_gyrobifastigium)) |
| Flattened gyrobifastigium (plesiohedron) | 8 | isosceles right triangles + silver rectangles | no (combinatorially F8-203) |
| Goldberg's octahedra | 8 | "at least 49" space-filling octahedra (Goldberg, Geom. Dedicata 10 (1981) 323–335); full list not accessible online | only 8-VI and the gyrobifastigium are known to us to have equilateral members; the full paper is worth checking |
| Goldberg's heptahedra | 7 | 16 of the 34 heptahedra fill space | only prisms are known to us to be equilateral (F7-031). F7-029 is triangle ⊕ triangle |
| Catoptric cells: tetrahedrille, pyramidille, oblate octahedrille (square bipyramid), … | 4–12 | fundamental domains of the cubic reflection groups ([Stereohedron](https://en.wikipedia.org/wiki/Stereohedron)) | only the cube and the RD. An equilateral square bipyramid is the regular octahedron, which does not tile |
| Ten-of-diamonds (Goldberg 10-II) | 10 | 8 isosceles triangles + 2 rhombi, combinatorially a square antiprism ([Wikipedia](https://en.wikipedia.org/wiki/Ten-of-diamonds_decahedron)) | equilateral members are rhombic or square antiprisms. In the symmetric ansatz they are rigid points; unlikely, but check (§9) |
| Half truncated octahedron (Goldberg 12-VIII); Goldberg 13-IV; 14-IV | 12–14 | surface-area candidates ([Ghang–Martin–Waruhiu](https://arxiv.org/abs/1305.1590)) | not equilateral (cut edges) |
| **Josehedron** (Fischer–Koch S TPMS extrema; Voronoi) | 12 | 4 isosceles triangles + 8 mirror-symmetric quadrilaterals, 12 vertices, 22 edges, 3 edge lengths ([Bernhard 2026](https://arxiv.org/abs/2604.07160)) | not equilateral. Its combinatorial type is a candidate (§9) |
| TPMS-derived cells (same paper, Table 1) | 6–20 | cube, TO, RD, triakis truncated tetrahedron, Laves 17-hedron, plus new 14-, 17- and 20-face types | only the classical ones are equilateral |
| Engel's stereohedra | 17–38 | 38 faces, Dirichlet stereohedron for a cubic group (Engel 1980–81) | — |
| Bounds | — | Dirichlet stereohedra: ≤ 25 facets for the full cubic groups and ≤ 92 for all 3D groups ([Sabariego–Santos III](https://arxiv.org/abs/math/0608039), [IV](https://arxiv.org/abs/0708.2114)). Wikipedia, citing Schmitt's 2016 thesis, gives 38 as the plesiohedron maximum | — |

## 6. Biology

- **Honeybee cell**: a hexagonal prism capped by 3 rhombi (or by Fejes Tóth's cap of 2 hexagons + 2 rhombi). Double combs stacked flat face to flat face tile space with one cell shape, 10 faces. The lateral faces are trapezoids with unequal parallel sides, and an equilateral version is impossible (§10). References: Fejes Tóth, *What the bees know and what they do not know*, Bull. AMS 70 (1964); D'Arcy Thompson, *On Growth and Form*.
- **Plant parenchyma and pith** (Matzke): irregular 14-hedra on average; not monohedral.
- **Epithelia**: prisms and "scutoids" (Gómez-Gálvez et al., Nat. Commun. 2018). Scutoids have curved faces and are not convex polyhedra.
- No biological source gave a new convex monohedral tile.

## 7. Cross-check against our ≤ 8-face list

| Natural or classical item | Our family |
|---|---|
| cube, rhombohedra (sc, pyrochlore Voronoi, golden rhombohedra) | F6-007 |
| triangular prisms (triangular prismatic honeycomb) | F5-001 |
| hexagonal prism (hexagonal WS cell, CHA *d6r*) | F8-192. Note: tilers include right prisms over **non-centrally-symmetric** equilateral Reinhardt type-2 hexagons (A + B + D = 360°). Example angles A…F = 93.724, 157.617, 103.242, 108.659, 138.100, 118.658; opposite angles differ by about 15–20°. By contrast A + B + C = 360° forces central symmetry (§10) |
| pentagonal prisms over plane-tiling pentagons | F7-031. Members whose base has two *adjacent* angles summing to 180° are rhombus ⊕ triangle with a rhombus side ∥ the triangle plane |
| gyrobifastigium, flattened gyrobifastigium, equilateral SCD biprism | F8-203 = triangle ⊕ triangle (8-face chamber) |
| boundary of F8-203 | F7-029 = triangle ⊕ triangle with one side ∥ the other triangle's plane: 3 T + 3 R + 1 pentagon, matching `5RRRTTT` |
| rhombus ⊕ triangle with a triangle side ∥ the rhombus plane | F8-249 |
| gabled rhombohedron / elongated gyrobifastigium (Goldberg 8-VI) | F8-255 |
| regular octahedron, square pyramid, triangular bipyramid, deltahedra | non-tilers (rigid; angle and Dehn filters) |

**Verdict: no equilateral convex ≤ 8-face space-filler from these sources is missing.**

## 8. Minkowski sums of unit polygons (observation from this survey)

**Dehn invariant.** For a Minkowski sum of polygons and segments in general position, every edge is a translate of an edge of one summand. For each such edge direction, the exterior dihedral angles sum to π. So Σθ over the edges parallel to that direction is a rational multiple of π, and D = 0 automatically. These sums therefore always pass the Dehn filter, and equal-length summand edges make them equilateral.

| Sum (unit edges) | Generic F | Face types | Dim | Tiles? |
|---|---|---|---|---|
| triangle ⊕ segment | 5 | prism | — | yes (F5-001) |
| triangle ⊕ triangle, chamber A | 8 | 4 T + 4 R; vertex degrees 4⁴3⁴ | 3 | yes: = F8-203 |
| **triangle ⊕ triangle, chamber B** | **9** | **4 T + 5 R**; V = 9 (5×deg 4, 4×deg 3), E = 16 | 3 | **yes, new**: lattice ⟨a₁+b₁, a₁+b₂, a₂−a₁⟩ plus the point reflection x ↦ 2(o₁+o₂) − x, for a suitable labelling T₁ = conv(o₁, o₁+a₁, o₁+a₂), T₂ = conv(o₂, o₂+b₁, o₂+b₂). Exact overlap check: 14/14 random members |
| rhombus ⊕ triangle, special positions | 7, 8 | pentagonal prism; F8-249 | — | yes |
| **rhombus ⊕ triangle, generic** | **10** | **2 T + 8 R**; V = 11, E = 19 | 4 | **yes, new**: lattice ⟨r₁+r₂+t₁, r₁+r₂+t₂, r₁−r₂+t₁−t₂⟩ plus x ↦ 2(o₁+o₂) + r₁ + t₁ − x, for a suitable labelling. Exact check: 12/12 members |
| rhombus ⊕ rhombus | 12 | 12 R | 5 | yes (RD-type parallelohedron) |
| triangle ⊕ centrally symmetric hexagon | 13 | 2 T + 2 hexagons + 9 R | 5 | not found (limited lattice + point-reflection search, 4 members) |
| triangle ⊕ triangle ⊕ segment | 14–15 | T + R | 5 | not found (limited search, 3 members) |

The "exact check" works as follows:
1. the volume equals |det L|/2;
2. no translate or reflected translate within ±3 lattice steps overlaps the tile in positive volume, which was checked with LP and half-space intersection.

These are numerical results and should be re-certified with the project's own `tilecheck.py` pipeline. As far as we found, the 9- and 10-face sums do not appear in the literature. Both are tilings by translates and point reflections, the class in Kuperberg's second question.

## 9. Candidates with 9–14 faces for the extension

| Priority | Candidate | F | Status |
|---|---|---|---|
| A | triangle ⊕ triangle (chamber B) | 9 | equilateral; tiles (§8) |
| A | rhombus ⊕ triangle (generic) | 10 | equilateral; tiles (§8) |
| A | RD-type zonotopes, 4 unit generators (fcc WS, Bilinski) | 12 | equilateral; tile |
| A | elongated dodecahedron, equilateral | 12 | equilateral; tiles |
| A | truncated octahedron family (bcc WS, Kelvin, sodalite) | 14 | equilateral; tiles |
| B | α-U Voronoi type (ED-like, not centrally symmetric) | 12 | natural tiler; equilateral realizability unknown |
| B | Josehedron type (4 triangles + 8 quadrilaterals, 12 vertices) | 12 | natural tiler; unknown |
| B | Williams β-tetrakaidecahedron (2 squares + 8 pentagons + 4 hexagons) | 14 | monohedral foam; unknown |
| B | the 10 simple 14-face isohedral tiles of Delgado-Friedrichs & O'Keeffe 2005 | 14 | combinatorial tilers; test each |
| B | other polygon sums: triangle ⊕ hexagon (13), triangle ⊕ triangle ⊕ segment (14), pentagon ⊕ triangle, … | 11–15 | D = 0 automatically; tiling open |
| C | ten-of-diamonds / rhombic antiprism | 10 | equilateral points are rigid; likely fail the angle filter |
| C | lonsdaleite / ice Ih cell | 11 | not equilateral; equilateral realizability unknown |
| — | trapezo-rhombic dodecahedron (hcp), bee cell | 12, 10 | **excluded** (§10) |
| — | triakis truncated tetrahedron (diamond) | 16 | equal edges collapse it to the regular tetrahedron |

## 10. Short proofs

**(a) Trapezo-rhombic dodecahedron: no convex equilateral realization.** The type consists of an upper cap of 3 quadrilaterals around a degree-3 apex, a belt of 6 quadrilaterals, and a lower cap of 3 quadrilaterals.
1. Equal edges make every face a rhombus.
2. Consecutive belt faces share their "vertical" edges, so all 6 vertical edges are equal to one vector w.
3. Hence the lower boundary hexagon is the upper one translated by −w.
4. In the hcp type, each vertical edge joins the *far* corner of an upper rhombus to the *far* corner of a lower rhombus. (In the rhombic dodecahedron it joins a far corner to a side corner.)
5. A cap of three rhombi is determined by its boundary: apex = s₁ + s₂ − f₁₂.
6. So the lower cap is the upper cap translated by −w. It bulges the same way, and the solid is not convex.

**(b) Bee cell** (flat hexagon, 6 lateral quadrilaterals, cap of 3 rhombi).
1. Equal edges turn the lateral faces into rhombi.
2. So all 6 lateral edges equal one vector w, and the cap's boundary hexagon is the base hexagon + w, which is planar.
3. Each cap rhombus has three of its vertices (side, far corner, side) on that boundary. They are not collinear, so the rhombus lies in the boundary plane. Its fourth vertex, the apex, must then also lie in that plane, and the cap is flat. Contradiction.

**(c) Equilateral hexagon with A + B + C = 360°.**
1. The turning angles at A, B, C sum to 180°, so CD = −FA.
2. Then AB + BC = −(DE + EF). Two unit vectors with the same nonzero sum are equal up to order.
3. The order AB = −EF, BC = −DE violates the monotonicity of edge directions, so AB = −DE and BC = −EF.
4. The hexagon is therefore centrally symmetric.

Type 2 (A + B + D) has no such constraint. Random solutions are non-centrally-symmetric (§7).

**(d) Triakis truncated tetrahedron.** With equilateral triangles, each cap over a truncation triangle is a regular tetrahedron. Its faces are coplanar with the adjacent hexagons, so the solid is the regular tetrahedron.

## References

- Wikipedia: [Wigner–Seitz cell](https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell), [Plesiohedron](https://en.wikipedia.org/wiki/Plesiohedron), [Stereohedron](https://en.wikipedia.org/wiki/Stereohedron), [Space-filling polyhedron](https://en.wikipedia.org/wiki/Space-filling_polyhedron), [Trapezo-rhombic dodecahedron](https://en.wikipedia.org/wiki/Trapezo-rhombic_dodecahedron), [Triakis truncated tetrahedron](https://en.wikipedia.org/wiki/Triakis_truncated_tetrahedron), [Elongated gyrobifastigium](https://en.wikipedia.org/wiki/Elongated_gyrobifastigium), [Ten-of-diamonds decahedron](https://en.wikipedia.org/wiki/Ten-of-diamonds_decahedron), [Einstein problem (SCD tile)](https://en.wikipedia.org/wiki/Einstein_problem), [A15 phases](https://en.wikipedia.org/wiki/A15_phases).
- [MathWorld, Space-Filling Polyhedron](https://mathworld.wolfram.com/Space-FillingPolyhedron.html): Goldberg's counts of 27 hexahedra, 16 of 34 heptahedra, and at least 49 octahedra.
- M. Goldberg: The space-filling pentahedra, JCTA 13 (1972) and II, JCTA 17 (1974); Three infinite families of tetrahedral space-fillers, JCTA 16 (1974); On the space-filling hexahedra, Geom. Dedicata 6 (1977); …heptahedra, 7 (1978); Convex polyhedral space-fillers of more than twelve faces, 8 (1979); On the dodecahedral space-fillers, 10 (1981); On the space-filling octahedra, 10 (1981) 323–335; On the space-filling enneahedra, 12 (1982); On the space-filling decahedra, Structural Topology (1982).
- O. Delgado-Friedrichs, M. O'Keeffe, [Isohedral simple tilings: binodal and by tiles with ≤ 16 faces](https://doi.org/10.1107/S0108767305009578), Acta Cryst. A61 (2005) 358–362.
- O. Delgado-Friedrichs, D. Huson, [4-regular vertex-transitive tilings of E³](https://link.springer.com/article/10.1007/s004540010035), DCG 24 (2000).
- O. Delgado-Friedrichs, M. O'Keeffe, D. Proserpio, M. Treacy, [Three-periodic nets, tilings and surfaces](https://doi.org/10.1107/S2053273323000414), Acta Cryst. A79 (2023) 192–202.
- N. Anurova, V. Blatov, G. Ilyushin, D. Proserpio, [Natural tilings for zeolite-type frameworks](https://doi.org/10.1021/jp1030027), J. Phys. Chem. C 114 (2010).
- [RCSR](https://rcsr.anu.edu.au); [EPINET](https://epinet.anu.edu.au).
- M. Baake, D. Frettlöh, [SCD patterns have singular diffraction](https://arxiv.org/abs/math-ph/0411052), J. Math. Phys. 46 (2005).
- R. Sabariego, F. Santos, On the number of facets of 3D Dirichlet stereohedra [III](https://arxiv.org/abs/math/0608039) and [IV](https://arxiv.org/abs/0708.2114).
- W. Ghang, Z. Martin, S. Waruhiu, [Surface-area-minimizing n-hedral tiles](https://arxiv.org/abs/1305.1590).
- M. Bernhard, [The Josehedron](https://arxiv.org/abs/2604.07160) (2026).
- B. Klaassen, [Which equilateral convex polygons tile the plane?](https://arxiv.org/abs/2506.18473), J. Geom. Graphics 29(2).
- A. Schoen, [Double Laves Voronoi cell](https://www.schoengeometry.com/e70-tpms-media/double_Laves_voronoy_cell.pdf).
- J. Lagarias, D. Moews, Polytopes that fill ℝⁿ and scissors congruence, DCG 13 (1995).
- L. Fejes Tóth, What the bees know and what they do not know, Bull. AMS 70 (1964) 468–481.
