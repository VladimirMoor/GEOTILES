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
