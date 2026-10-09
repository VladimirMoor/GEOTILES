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
