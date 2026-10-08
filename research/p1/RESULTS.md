# Problem 1: equilateral convex polyhedra with ≤ 8 faces that tile space

First pass, 2026-10-07. The class covers every convex polyhedron with all edges of length 1 and 4–8 faces, and allows any monohedral tiling: copies may use any isometries, and the tiling need not be face-to-face or isohedral.

## Pipeline

| Step | Script | Result |
|---|---|---|
| 1. Combinatorial types | `enumerate_types.py` (plantri `-pc3 -d`) | 1 + 2 + 7 + 34 + 257 = **301** types (OEIS A000944) |
| 2. Equilateral convex realizations | `realize.py` | **35** types are realizable: 16 rigid types (17 solids) and 19 continuous families |
| 3. Rigid solids | `analyze.py` | **none tile** |
| 4. Families: Dehn invariant | analysis below | **11** families are excluded entirely |
| 5. Families with zero Dehn invariant | `identities.py`, `search_tiling.py`, `tilecheck.py` | **5** families where every member tiles; prisms are partly open |

### Step 2 in detail

The unknowns are the vertex positions, the face planes and 6 gauge conditions. The equations are |x_u − x_v|² = 1 on every edge, n_f·x_v = d_f for every vertex of a face, and |n_f| = 1. By Euler's formula this system is square, so a generic type has finitely many solutions. Flexible types instead form continuous families, whose dimension we read off as the nullity of the Jacobian minus 6.

Solutions come from Levenberg–Marquardt runs with 2000 starts per type. Each solution is verified strictly:
- edge lengths equal 1 to within 1e-7;
- every face plane is supporting, with all other vertices at least 1e-6 below it;
- every face is a strictly convex polygon in the cyclic order of the combinatorial type.

**Non-degeneracy.** A type counts only if some solution has every dihedral angle at least 0.5° away from 0° and 180°. Genuine types have a margin of at least 4.4°. Two types, F7-011 and F8-227, produced only near-boundary artifacts with margins below 0.001° (two faces practically coplanar). After Newton refinement these artifacts are not strictly convex, so both types are excluded.

**Independent rerun.** A second run with 6000 starts per type and a different seed found exactly the same 35 types, the same 17 rigid solids with identical angles (each hit 47–2916 times), and the same flexibility dimensions.

**Caveat.** This is still a numerical search, not a proof that the list is complete.

## Results

### Rigid solids (17): none tile

- **15 fail the dihedral angle filter.** Around an edge of a tiling the angles must sum to 360°, or to 180° when the edge lies on a neighbour's face. These are the regular tetrahedron, J1, J2, J12, the octahedron and 10 non-regular-faced solids.
- **Two pass the filter but still do not tile:** the truncated tetrahedron (F8-199) and the octahedron + tetrahedron (F7-030, which is the 60° rhombohedron minus one corner tetrahedron). Both decompose into unit regular tetrahedra and octahedra.
  - In any tiling by such pieces, the angles around an edge satisfy a·arccos(1/3) + b·(π − arccos(1/3)) (+π) = 2π. Since arccos(1/3)/π is irrational, this forces a = b.
  - Counting edge length then forces N_tet = 2 N_oct.
  - F7-030 has the ratio 1 : 1 and the truncated tetrahedron has 7 : 4, so neither can tile.

### Dehn invariant lemma

For an equilateral polyhedron, D(P) = Σ ℓ_e ⊗ θ_e = 1 ⊗ Σθ_e, which vanishes iff **Σθ_e ∈ ℚπ**. By Lagarias–Moews (1995), any polytope that tiles ℝ³ with congruent copies has D = 0.

Along every family Σθ is constant (spread < 1e-11). Integer relation detection (PSLQ) identifies each sum as follows:

| Family | Faces | Dim | Σθ | Verdict |
|---|---|---|---|---|
| F6-002 | 5RRTTT | 1 | 3π + 6·arccos(1/3) | does not tile |
| F7-009 | RRRTTTT | 2 | 4π + 6·arccos(1/3) | does not tile |
| F7-016 | 6RRTTTT | 1 | 9π − 6·arccos(1/3) | does not tile |
| F8-057 | 5RRTTTTT | 1 | 3π + 12·arccos(1/3) | does not tile |
| F8-178 | 7RRTTTTT | 1 | 3π + 5·arccos(1/√5) + 15·θ(J2) | does not tile* |
| F8-204 | 6RRRTTTT | 1 | 6π + 6·arccos(1/3) | does not tile |
| F8-206 | RRTTTTTT | 1 | 10π − 6·arccos(1/3) | does not tile |
| F8-235 | 5RRRRTTT | 1 | 11π − 6·arccos(1/3) | does not tile |
| F8-244 | 556RRRTT | 1 | 7π + 6·arccos(1/3) | does not tile |
| F8-250 | 5RRRRTTT | 1 | 11π − 6·arccos(1/3) | does not tile |
| F8-252 | RRRRTTTT | 2 | 6π + 6·arccos(1/3) | does not tile |

Here R = rhombus, T = equilateral triangle, and a digit is an n-gon. These families are "regular-faced cap + oblique prism". Prisms have Dehn invariant 0, so D(family) = D(cap) ≠ 0. *The irrationality for the J2 cap still needs a written proof.

### Families with D = 0

| Family | Dim | Description | Tiles? |
|---|---|---|---|
| F5-001 | 2 | triangular prisms (all slants) | every member: half of a parallelepiped |
| F6-007 | 3 | rhombic parallelepipeds | every member: lattice tiling |
| F8-203 (+ F7-029 on its boundary) | 3 | **generalized gyrobifastigium**: two oblique triangular prisms on a common rhombus | every member: lattice a, b, apex₁ − apex₂ + point reflection x ↦ a + b + apex₁ − x |
| F8-249 | 3 | **rhombus ⊕ triangle** (Minkowski sum; one side of the triangle parallel to the rhombus) | every member: lattice d₁, a + b + d₂, b − a + point reflection in the midpoint of a rhombus edge |
| F8-255 | 3 | **generalized elongated gyrobifastigium**: rhombic parallelepiped + two triangular prisms whose triangles merge with its faces into pentagons | every member: lattice a, b, c + apex₁ − apex₂ + point reflection |
| F7-031 | 4 | equilateral pentagonal prisms | **open**; tiles if the cross-section ⊥ the lateral edge tiles the plane |
| F8-192 | 5 | equilateral hexagonal prisms | **open**; same sufficient condition |

The families marked "every member" are backed by explicit tilings checked numerically. Volume ratio V(tile)/V(fundamental domain) = 1, and Monte Carlo coverage passed for 20 random members of F8-203, 18 of F8-255 and 40 of F8-249 in the browser, plus 30 more in Python. For every tile-able family, every edge class is covered by an identity of the angle filter that holds over the whole family.

**Prism remark (sufficient condition).** An oblique prism tiles space whenever its cross-section perpendicular to the lateral edges tiles the plane. Cut the infinite rods of any plane tiling by that cross-section; each rod is a congruent copy of the original, so each can be cut into congruent prisms independently. The cross-section is an affine image of the equilateral base, and every affine image (up to similarity) occurs for some slant.

## Open items

1. Completeness of step 2: confirmed by an independent 6000-start rerun; a certified method (homotopy continuation or interval arithmetic) is still needed for a proof.
2. A rigorous structural identification of each family, replacing the numerical one.
3. Pentagonal and hexagonal prisms: which members tile. Is the sufficient condition also necessary?
4. Check the literature for the three "generalized" families.

## Prisms (second pass, 2026-10-08)

**Context.** Whether the base of a right prism that tiles space must itself tile the plane is an open problem. W. Kuperberg posed it at the 1993 Smith College tiling session; it appears in Senechal's problem list and on [Eppstein's page](https://ics.uci.edu/~eppstein/junkyard/tiling-problems.tex). Settling equilateral prisms completely would settle special cases of this question.

Kuperberg's second question from the same session asks for the convex polyhedra that tile with translates of themselves and of their point reflection. All five of our tile-able families tile in exactly this way.

### Right equilateral pentagonal prisms (`prisms.py`)

The dihedral angles are the pentagon angles γᵢ and 90°. The filter requires, for each i, Σ aⱼγⱼ ∈ {360°, 270°, 180°, 90°} with aᵢ ≥ 1.

| Stratum | Dim | Base tiles the plane? | Prism tiles space? |
|---|---|---|---|
| two angles sum to 180° | 1 | yes (Klaassen 2025: any two angles) | **yes** (stack layers) |
| two **non-adjacent** angles sum to 270° | 1 | no | **open** (passes the filter) |
| two adjacent angles sum to 270° | — | not realizable by an equilateral convex pentagon | — |
| isolated points (relations with Σaⱼ ≤ 7) | 0 | only P7 tiles; P7 was rediscovered automatically | P7: yes; **5 others open** |

The five isolated open pentagons have angles (degrees, in cyclic order A–E):
- (141.332, 77.337, 126.847, 106.307, 88.178)
- (91.580, 73.982, 176.841, 44.210, 153.387)
- (50.436, 166.805, 79.129, 92.324, 151.307)
- (34.418, 177.061, 76.745, 79.684, 172.092)
- (36.000, 173.010, 78.990, 78.990, 173.010)

The point search is complete only for relations with Σaⱼ ≤ 7.

### Hexagonal prisms: a simplification

**Lemma (numerical, to be proved).** An equilateral convex hexagon with alternating angle sum A + C + E = 360° is centrally symmetric. In 300 random samples, opposite angles were equal to 1e-6.

So the codimension-1 filter strata of right equilateral hexagonal prisms are:
- adjacent triples summing to 360°, which tile the plane (Klaassen);
- alternating triples, which give centrally symmetric hexagons. These are parallelohedra and tile for every slant.

Any residual is therefore at most 1-dimensional. Computing it is still to do.

### Oblique prisms

The lateral dihedral angles are the angles of the cross-section Q perpendicular to the lateral edge. Q is an affine image of the base, and Q is no longer equilateral, so Klaassen's criterion does not apply. Consider the filter stratum "two non-adjacent angles of Q sum to 180°". There Q passes the planar angle condition but in general fails the side-length condition of type 2. This stratum is 3-dimensional, so the open part of oblique prisms is large.

## Rigorous verification (in progress, 2026-10-08)

**Exact model.** The system is posed over ℚ (`exact.py`):
- unit edges;
- rhombi as parallelograms (linear equations);
- face planes through unit normals, |n|² = 1;
- Rabinowitsch variables for non-degenerate faces and for non-coplanar adjacent faces.

msolve computes Gröbner bases and isolates every real solution with rational intervals. Strict convexity and the dihedral filter are then decided with interval arithmetic (`rigid_exact.py`).

Each certified claim is a proof, modulo the correctness of msolve. Two pitfalls were found and fixed:
- msolve may permute variables, so the order is read from its `-P 1` output;
- its parser mishandles unexpanded expressions such as `-(-35)`, so every polynomial is expanded first.

Every solution read back from msolve is also checked against the equations in 260-bit arithmetic.

| Outcome over ℂ | Types | Status |
|---|---|---|
| no solutions at all | 135 | **proved: no convex realization** |
| finitely many solutions | 35 | **certified.** 28 have no convex real solution. 7 have exactly one convex realization: tetrahedron, J1, J12, J2, F7-021, F7-030, octahedron. Six of these fail the filter with interval proof; F7-030 passes it and is excluded by the tetrahedron–octahedron counting lemma |
| positive-dimensional | 48 | 14 are the known families. The other 34 are being decided by `curve1d.py`, a rigorous sweep of 1-dimensional components |
| Gröbner basis not finished in 300 s | 83 | rerunning with a 2 h limit (14 of these are numerically realizable) |

So far **170 of the 301 types are settled rigorously**, and no rigorous result contradicts the numerical enumeration.

### Final status of the rigorous pass (2026-10-08)

Machine-readable lists are in `data/rigor_status.json`.

**Certified (170 of 301 types).**
- 135 types have no complex solution at all, so they have no convex realization.
- 35 types have finitely many solutions, all certified. Seven of them have exactly one convex realization, and none of those seven tiles space.

**Not certified (131 types).**
- 48 have positive-dimensional components. 14 of these are the known families.
- 83 types: msolve did not finish within 300 s, and not even within 2 h for the first ones tried.

For these 131 types we rely on two independent numerical searches (2000 and 6000 starts, different seeds), which agree exactly. 103 of them have no numerical realization. The remaining 28 are the numerically realizable rigid solids and families, which are analysed as described above.

**Approaches tried for the hard types that did not work:**
- saturation in Singular: too slow;
- local "patch" subsystems: no contradiction is local;
- an edge-class (zone) formulation: it does not reduce the size, because the hard types have few rhombi;
- longer msolve runs: still unfinished after 2 h.

A complete proof would need a certified method suited to these sizes. Candidates are certified homotopy continuation (e.g. alpha-certified monodromy with a trace test) or interval branch-and-bound on a reduced parametrization.

The rigorous **1-dimensional sweep** (`curve1d.py`) works correctly. On F6-002 it finds the convex family members, but it is too slow for routine use: about 50 minutes per sample fibre.

## Minkowski sums of unit simplices (2026-10-08)

The literature survey (`NATURE_SURVEY.md`) noticed that our tiling families are Minkowski sums of unit polygons and segments, and reported two new equilateral tilers with 9 and 10 faces. We re-checked everything with our own code (`minkowski.py`, `search_tiling.py`, and the site's tiler).

### Confirmed equilateral space-fillers

All are numerical results: volume ratio 1 and Monte Carlo coverage, on dozens of random members each.

**Triangle ⊕ triangle, general position.** Two chambers, about 50/50 for random positions.
- 8 faces: the generalized gyrobifastigium (F8-203).
- **9 faces: 4 triangles + 5 rhombi, 9 vertices, 3 parameters. New.**

Writing T₁ = o₁ + {0, u₁, u₂} and T₂ = o₂ + {0, v₁, v₂}, the tilings are:

| Chamber | Lattice | Point reflection x ↦ w − x |
|---|---|---|
| 8 faces (formula A) | ⟨v₁, u₁, u₂ − v₂⟩ | w = 2(o₁+o₂) + v₁ + u₁ + u₂ |
| 9 faces (formula B) | ⟨u₁+v₁, u₁+v₂, u₂−u₁⟩ | w = 2(o₁+o₂) |

Each formula works for exactly 4 of the 36 labellings of a member, and fails on the other chamber.

**Rhombus ⊕ triangle, general position: 10 faces (2 triangles + 8 rhombi), 11 vertices, 4 parameters. New.**
- Lattice ⟨r₁+r₂+t₁, r₁+r₂+t₂, r₁−r₂+t₁−t₂⟩.
- Point reflection x ↦ 2(o+o_T) + r₁ + t₁ − x.
- Works for 4 of the 48 labellings.
- Its special position, with a triangle side parallel to the rhombus, is F8-249.

### Survey of sums: total simplex dimension decides

| Sum (general position) | Σ dim | F | Angle filter on random members | Tiling found |
|---|---|---|---|---|
| 3 segments, triangle ⊕ segment | 3 | 6, 5 | passes | yes |
| triangle ⊕ triangle | 4 | 8 / 9 | 6/6 pass | yes |
| rhombus ⊕ triangle (= 2 segments ⊕ triangle) | 4 | 10 | 6/6 pass | yes |
| rhombus ⊕ rhombus (4 segments) | 4 | 12 | passes | yes (translations) |
| centrally symmetric hexagon ⊕ segment | 4 | 8 | passes | yes (translations) |
| triangle ⊕ triangle ⊕ segment | 5 | 14–15 | **0/6** | no |
| triangle ⊕ hexagon | 5 | 13 | **0/6** | no |
| rhombus ⊕ triangle ⊕ segment | 5 | 17 | **0/6** | no |
| pentagon ⊕ triangle, pentagon ⊕ segment | — | 11–12, 7 | **0/6** | no |

A failed angle filter is a proof that the member does not tile in any way, so the negative entries are rigorous for the sampled members.

**Conjecture.** A Minkowski sum of unit segments and unit equilateral triangles in general position in ℝ³ tiles space if and only if the total dimension of the summands is at most 4 = d + 1. The tiling then uses only translations and point reflections.

This is analogous to McMullen's zonotope criterion: 4 generic generators tile, 5 do not. The planar analogue holds: triangle ⊕ segment is a pentagon with two parallel sides, which tiles.

**Question.** Is every equilateral convex space-filler a Minkowski sum of unit polygons and segments? It is true for every tiler we know with ≤ 8 faces: prisms, parallelohedra, the generalized gyrobifastigia and F8-249. It is also true for all the parallelohedra.

A natural route to a proof of the conjecture is the Cayley trick: sums of simplices are projections of products of simplices, and their subdivisions correspond to triangulations of the Cayley polytope.
