# Problem 2: new space-filling tetrahedra

Started 2026-10-09.

## Literature check

**No complete classification of space-filling tetrahedra exists.**

- **Classical tilers.**
  - Sommerville (1923) found four tetrahedra, No. 1–4; they tile face-to-face, orientation-preserving. Edmonds (2007) showed the list of such tiles is complete.
  - Hill (1896) and Baumgartner (1971) found further examples.
  - Goldberg (1974) found three one-parameter families. His first family is Hill's family 𝓕₁.
- **Recent work.**
  - Bongiovanni–Diaz–Kakkar–Sothanaphan, arXiv:1709.04139: sorted all tetrahedra into 25 edge-length types and identified the tiles among those with all dihedral angles of the form 2π/n.
- **Rational tetrahedra** (all dihedral angles in ℚπ, so Dehn invariant 0). Kedlaya–Kolpakov–Poonen–Rubinstein (2020) classified them: two one-parameter families 𝓕₁ (Hill) and 𝓕₂, plus 59 sporadic tetrahedra.
- **Chentouf–Sun, "Tetrahedra tiling problem", arXiv:2312.01654 (Dec 2023).**
  - 𝓕₁ tiles; exactly one member of 𝓕₂ tiles.
  - At most 40 sporadic tetrahedra tile. These form the set 𝒜, and the paper leaves them **undecided**.
  - This disproves the converse of Debrunner's theorem.
- **Later work.** Searching on 2026-10-09 found no work that settles 𝒜. The companion paper arXiv:2312.01282 studies Dehn-invariant-zero tetrahedra in general.

**Target chosen: decide the 40 sporadic tetrahedra of 𝒜.** Any member shown to tile would be a new space-filling tetrahedron.

Data: `tetra.py`. Angles are listed in the order (α12, α34, α13, α24, α14, α23), as multiples of π.

## Results so far

### 1. The printed list 𝒜 has 42 entries, not 40

Two of the printed tuples are the ones the authors themselves exclude by the LP criterion in §3.1:
- (3/20, 11/20, 11/20, 1/4, 1/3, 2/3), our #34;
- (11/60, 31/60, 11/20, 1/4, 3/10, 7/10), our #35.

The true 𝒜 is the other 40.

All 42 tuples are verified to be genuine tetrahedra:
- the Gram matrix is singular to within 1e-40, PSD, with a positive kernel;
- the dihedral angles recomputed from the vertices agree to within 1e-60.

### 2. Their edge-LP reproduced exactly (`obstruct.py`)

Corollary 2.13 excludes exactly #34 and #35, and passes all 40 members of 𝒜. This independently checks both the transcription and the implementation.

### 3. Two members of 𝒜 are classical tilers (`known.py`, `isohedral.py`)

| # | Angles | Identification |
|---|---|---|
| 1 | (1/4, 1/3, 1/3, 1/4, 1/2, 2/3) | **Sommerville No. 3** (Table 1 of arXiv:1709.04139, after relabelling) |
| 3 | (1/4, 1/2, 1/2, 1/3, 1/3, 1/2) | **Sommerville No. 2** = Coxeter simplex [4,3^{1,1}] (all four faces are mirrors) |

The isohedral search (`isohedral.py`) checks the Poincaré polyhedron theorem over every face-pairing scheme. It finds 6 tilings for #1 and 4 for #3.

No other member of 𝒜 is a rational member of Goldberg's second or third family. Those families have only three rational members, and none of them is sporadic.

**So the genuinely open set has 38 tetrahedra.**

### 4. None of the 38 tiles face-to-face (`stars.py`)

**Edge-star criterion.** In a face-to-face tiling, the tiles around an edge segment PQ form a closed cycle of wedges. Consecutive wedges share the same triangle on base PQ, and the angles sum to 2π.

Triangles on base PQ are described by the pair of distances (|P·apex|, |Q·apex|). The existence of such a cycle through T's own edge e is then a closed-walk problem in a small graph. Mirror copies are included automatically.

**Result.** For each of the 38 open tetrahedra, some edge has no star. Hence **none of them admits a face-to-face tiling**. Control: #1 and #3 have stars on every edge.

**Rigour.**
- Angles are exact rationals.
- Length classes are separated by at least 0.0093 and coincide to within 1e-59.
- Wrongly merging two lengths would only add transitions, so non-existence is safe in that direction.

**What is new compared with the paper.** Chentouf–Sun's Proposition 2.7 already rules out face-to-face tilings for 26 of the 38. It is new for **12**: #2, 5, 7, 12, 19, 20, 22, 23, 25, 27, 36, 42.

| Open tetrahedra | Count |
|---|---|
| #2, 4–33 except 34, 35, and 36–42 | 38 |
| of which face-to-face exclusion is new | 12 |

### 5. What remains: tilings that are not face-to-face

Every one of the 40 has π-combinations of dihedral angles, so Proposition 2.3 never excludes non-face-to-face tilings (`nonf2f.py`).

**Structure found** (`halfstars.py`). A "bad" edge has no metric 2π-star. Bad edges also have no metric π-half-star. Consequently, at every point of every copy of a bad edge, in any tiling, the cycle of tiles has at least two "faults":
- either a half-plane in which the faces of the two adjacent tiles differ,
- or the edge lies in the interior of another tile's face.

In particular, every tiling by one of the 38 tetrahedra is non-face-to-face in an essential and dense way: faults occur along every copy of every bad edge.

**Weak criteria.**
- Vertex LP (`obstruct.py`): solid angles are rational, but the LP is too weak; all pass.
- Face-angle relations (`faceangles.py`): the face angles admit many relations summing to π or 2π, so planar vertex conditions alone are weak.

## Next steps (ideas)

1. **Fault-plane analysis.** In a non-face-to-face tiling, the faces lying in a fault plane Π form two different planar tilings of the same region S, one from each side. Along every interior edge of each planar tiling there is a π-combination containing the two adjacent tiles' angles, and interior vertices satisfy planar angle conditions. Develop this into a finite search.
2. **Constructive side.** Search for non-face-to-face tilings, for example dissections of rods or slabs whose cross-section tiles the plane. Since faults are forced along bad edges, the bad edges would have to lie on the sliding planes.
3. **Exact certification** of the length classes in sympy, as minimal polynomials, for a publication-grade write-up.

## Fault-plane analysis (2026-10-09, second pass)

**Planar conditions** (`faults.py`). Assume a tiling that is not face-to-face, and let Π be a fault plane. The faces lying in Π give two different planar tilings 𝒫₊ and 𝒫₋ of the same region S, one from each side. Two necessary conditions follow:

- **(E)** Along an interior edge of 𝒫₊ with adjacent triangle edges g and g′, the tiles above form a half-star of total π. Hence π − α_g − α_{g′} must be a non-negative integer combination of dihedral angles.
- **(V)** At an interior vertex, the face angles sum to 2π. At a T-junction (the vertex lies on another triangle's edge) they sum to π on one side.

**Result: these local conditions do not exclude anything.** All 38 open tetrahedra have both full and T-junction planar vertex stars. Example: six copies of one face around a point, which is the standard triangle tiling of the plane.

**Structural lemmas:**
- **Bad edges always lie in fault planes.** At every point of every copy of a bad edge (no metric 2π-star, no metric π-half-star), the cycle of tiles splits into metric chains separated by mismatched half-planes, or there is a face-interior event. Each mismatched half-plane lies in a fault region of its plane.
- **Bad edges sit on cluster edges.** Suppose a tiling comes from a convex cluster Q of copies glued face-to-face inside, and Q itself tiles. Then interior edges of Q need full metric stars, and edges lying in faces of Q need metric half-stars. So **every copy of a bad edge lies on an edge of Q**, and Q's dihedral angle there is the sum of a metric chain through that edge. The admissible chain sums are listed by `chains.py`. They are very restricted: for #17, #19, #25 a bad edge's chain is just the tile itself.

**Constructive search** (`clusters.py`, `poincare_q.py`, `run_clusters.py`).
- **Method.** We enumerate all clusters of k ≤ 7 copies glued face-to-face (mirror copies allowed), up to congruence, and keep the convex ones. For each convex cluster we run the Poincaré theorem with whole polygonal faces paired by isometries: translations, rotations, reflections and point reflections.
- **Validation.** For #1 and #3 there are 31 and 28 convex clusters, and many of them tile.
- **The 38 open tetrahedra.** They have only 0–24 convex clusters, a count that hardly grows with k, and **none of these clusters tiles**. Logs: `data_clusters5.log`, `data_clusters7.log`.

**Status.** The 38 open tetrahedra are still undecided. Any tiling would have to be non-face-to-face, with faults along every copy of every bad edge. It also could not come from a convex cluster of at most 7 copies tiling isohedrally.

## Attempt at #19 and #25 (2026-10-09)

#19 = (1/5, 1/5, 2/5, 1/3, 1/2, 2/3) and #25 = (1/5, 1/3, 1/3, 1/5, 1/2, 4/5). In both, the only bad edge is 23. It has a unique length and the largest dihedral angle (2π/3 and 4π/5).

**Forced combinations at the bad edge 23** (exact enumeration):
- **#19.** At every point of every copy of edge 23, the cycle is either {23, 23, 23}, or it contains a wedge at edge 24 (angle π/3). Every combination with exactly one 23 contains 24.
- **#25.** At every point there is a wedge of angle π/5 at edge 12 or 24, whose length is φ·|23|. The possible cycles are:
  - 23 + 6×(π/5);
  - 23 + (π/5) + 3×(π/3);
  - 23 + 23 + 2×(π/5);
  - 23 + (π/5) + face.
- **Faults.** At each point of the edge, at least one of the two faces of the tile at 23 is mismatched.

**Why local methods cannot succeed here.**
- **Edge LP.** It is feasible, and no combination is forced: every combination has minimum weight 0 (`lpcore.py`).
- **Icosahedral structure.** Both tetrahedra are icosahedral. All dihedral angles lie in {π/5, π/3, π/2, 2π/5, 2π/3, 4π/5}, and every vertex figure is a union of 1–14 Möbius triangles (2,3,5) of the icosahedral group H₃. Vertex 4 of #19 and vertices 1 and 4 of #25 *are* the (2,3,5) triangle, which tiles the sphere by reflections.
- **Consequence.** Spherical vertex stars always exist, and planar fault stars exist as well (see above). Every local angular condition is compatible with the H₃ mirror arrangement.

**The obstruction, if there is one, is metric and global.**
- The edge lengths involve the golden ratio:
  - #19: |24|/|23| = φ;
  - #25: |12|/|23| = |24|/|23| = φ.
- The distinct lengths are numerically linearly independent over ℚ: PSLQ finds no relation with coefficients up to 10⁶.

**Proposed global route (not yet carried out):**
1. A 1-D length balance along maximal segments of lines that bound tile faces. The endpoints of such a segment are tile vertices on both sides. Because the lengths are ℚ-independent, this forces equal counts of each edge type on the two sides, apart from face-interior parts.
2. Combine this with the forced wedges at edge 23 (for #25, a φ-length edge always accompanies a unit-length edge), aiming for a contradiction by counting.
3. Alternatively, show that the tiles use finitely many orientations (the H₃ orbit), unless some fault plane is a full plane. Then apply translation-type (Hadwiger) invariants.

## Dissection relations and the golden box (2026-10-09)

`halves.py`, `dissections.py`.

**Halves.** #25, #27 and #42 have a mirror symmetry (1 ↔ 4). In each case the mirror cuts the tetrahedron into two congruent halves, and every half is **#36** = (1/5, 1/2, 1/2, 1/3, 2/5, 1/2).

Each of #25, #27, #42 is therefore made of two mirror copies of #36:
- if #36 does not tile, then none of #25, #27, #42 tiles;
- if any of #25, #27, #42 tiles, then #36 tiles.

The halves themselves are not excluded by the edge LP, by Proposition 2.3, or by edge stars.

Searching all tetrahedral clusters of up to 6 copies finds no other relations inside 𝒜. The controls behave as expected: two copies of #1 make Hill(1/4), and #3 relates to #1, Sommerville No. 1 and Hill(1/4).

**#36 is the corner of a golden box.**
- The three dihedral angles at vertex 3 are π/2, so the edges 31, 32, 34 are mutually perpendicular.
- Their lengths are in the ratio 1 : φ : φ², with φ the golden ratio.
- All four faces lie in mirror planes of the icosahedral group H₃. The fourth normal is ½(φ, 1, φ⁻¹).

**Golden box = 4 × #36 + 1 × #7.**
- The box [0,1]×[0,φ]×[0,φ²] splits into four corner tetrahedra at alternate vertices, all congruent to #36 up to mirror image.
- The remaining central disphenoid has angles (1/5, 1/5, 1/3, 1/3, 3/5, 3/5), which is exactly **#7**, with twice the volume of #36. Verified to 50 digits.
- **So #36 and #7 tile space together** (two prototiles, a periodic tiling). Whether either one tiles alone is open.
- #7 is not a union of copies of #36: its four faces are acute, while two copies of #36 always expose right triangles.
