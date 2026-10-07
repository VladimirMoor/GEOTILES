# Problem 1: equilateral convex polyhedra with ≤ 8 faces that tile space

First pass, 2026-10-07. The class covers every convex polyhedron with all edges of length 1 and 4–8 faces, and allows any monohedral tiling: copies may use any isometries, and the tiling need not be face-to-face or isohedral.

## Pipeline

| Step | Script | Result |
|---|---|---|
| 1. Combinatorial types | `enumerate_types.py` (plantri `-pc3 -d`) | 1 + 2 + 7 + 34 + 257 = **301** types (OEIS A000944) |
| 2. Equilateral convex realizations | `realize.py` | **36** types are realizable: 16 rigid types (17 solids) and 20 continuous families |
| 3. Rigid solids | `analyze.py` | **none tile** |
| 4. Families: Dehn invariant | analysis below | **11** families are excluded entirely |
| 5. Families with zero Dehn invariant | `identities.py`, `search_tiling.py`, `tilecheck.py` | **5** families where every member tiles; prisms are partly open |

### Step 2 in detail

The unknowns are the vertex positions, the face planes and 6 gauge conditions. The equations are |x_u − x_v|² = 1 on every edge, n_f·x_v = d_f for every vertex of a face, and |n_f| = 1. By Euler's formula this system is square, so a generic type has finitely many solutions. Flexible types instead form continuous families, whose dimension we read off as the nullity of the Jacobian minus 6.

Solutions come from Levenberg–Marquardt runs with 2000 starts per type. Each solution is verified strictly:
- edge lengths equal 1 to within 1e-7;
- every face plane is supporting, with all other vertices at least 1e-6 below it;
- every face is a strictly convex polygon in the cyclic order of the combinatorial type.

**Caveat.** This is a numerical search, not a proof that the list is complete. An independent rerun with 6000 starts and a new seed is in progress.

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
| F7-011 | — | only near-degenerate numerical solutions | to be checked |

The families marked "every member" are backed by explicit tilings checked numerically. Volume ratio V(tile)/V(fundamental domain) = 1, and Monte Carlo coverage passed for 20 random members of F8-203, 18 of F8-255 and 40 of F8-249 in the browser, plus 30 more in Python. For every tile-able family, every edge class is covered by an identity of the angle filter that holds over the whole family.

**Prism remark (sufficient condition).** An oblique prism tiles space whenever its cross-section perpendicular to the lateral edges tiles the plane. Cut the infinite rods of any plane tiling by that cross-section; each rod is a congruent copy of the original, so each can be cut into congruent prisms independently. The cross-section is an affine image of the equilateral base, and every affine image (up to similarity) occurs for some slant.

## Open items

1. Completeness of step 2: the 6000-start rerun, then possibly certified methods (homotopy continuation or interval arithmetic).
2. A rigorous structural identification of each family, replacing the numerical one.
3. Pentagonal and hexagonal prisms: which members tile. Is the sufficient condition also necessary?
4. F7-011.
5. Check the literature for the three "generalized" families.
