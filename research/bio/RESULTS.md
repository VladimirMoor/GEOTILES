# Protein crystals: Dirichlet cells of molecules (side project)

Question: in real protein crystals, how many neighbours does one molecule have? We measure this by the face count of the Dirichlet–Voronoi cell of the molecule's centre under the full space group. The record stereohedra of problem 4 (Engel 38, in I4₁32) are DV cells of exactly this kind, with a hand-picked point position. Here the positions come from nature.

Method (`pdb_dv.py`):
- Get cell, space group and Cα atoms of model 1 from RCSB (GraphQL plus ModelServer).
- Take the centre as the centroid of the deposited asymmetric unit.
- Convert to fractional coordinates and build the DV cell of its orbit with `research/p4/dvfast.py`. This is a float computation, not certified.

## I4₁32 (the group of Engel's 38-hedron)

- 156 X-ray entries, 146 analysed.
- 10 entries failed. They have no Cα atoms, i.e. they are nucleic acid or ligand-only entries.

| Faces | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|
| Entries | 9 | 31 | 86 | 10 | 8 | 2 |

The maximum is 20 faces (1JKY, 4PPI). The next are 19 faces: 1A87, 3LJE, 5A9W/X/Y and 5OID.

Real molecules sit at generic positions far from Engel's special point, and their cells are moderate: 17 faces is typical. Nature does not come near the 38 record. The centroid of the asymmetric unit is a crude proxy for "the molecule", since one ASU can hold several chains.

Run: `nice -n 10 ../../.venv/bin/python pdb_dv.py "I 41 3 2"`. The output is `data/dv_I4132.json`, and the download cache is git-ignored.

## I4₁22 (the group of Schmitt's 35-hedron)

- 1129 entries, 1117 analysed.
- The 12 most-faceted cells were re-certified in exact rational arithmetic with `research/p4/exact.py`. All 12 counts agree with the float computation.

| Faces | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 28 | 29 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Entries | 111 | 219 | 189 | 70 | 79 | 158 | 81 | 149 | 14 | 31 | 5 | 3 | 3 | 2 | 3 |

The maximum is **29** (7ZC0 at c/a 0.275; 8S97 and 8S9D at c/a ≈ 1.55). Next come 28 (2X2B, 5N7L) and 26 (4LMR, 7VOB, 7VOC). The record for the group is 35 (Schmitt).

## P6₁22

- 3328 entries, 3235 analysed.
- The top 12 were certified exactly, and all agree.

| Faces | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 24 | 30 | 31 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Entries | 92 | 178 | 318 | 868 | 970 | 406 | 223 | 84 | 56 | 36 | 1 | 1 | 2 |

The maximum is **31** (4UX6 and 8GS1, c/a ≈ 0.55), followed by 30 (6P64, c/a 0.41). The record for the group is 34 (Schmitt; our certified maximum is 33).

## Summary

| Group | Entries | Typical | Max in proteins | Record stereohedron |
|---|---|---|---|---|
| I4₁32 | 146 | 17 | 20 | 38 (Engel) |
| I4₁22 | 1117 | 15–21 | 29 | 35 (Schmitt) |
| P6₁22 | 3235 | 16–17 | 31 | 34 (Schmitt) |

In the non-cubic groups, a few natural crystals reach 83–91% of the record face count, and always in strongly anisotropic lattices. That is the same regime where the record stereohedra live. Most molecules have 15–18 contacting neighbours.

## Common space groups and single chains (2026-10-25)

**Questions.**
1. How many neighbours does a molecule have in the most common protein space groups, which have fewer symmetry operations than the groups above?
2. How does the count change when each protein chain gets its own point instead of the centroid of the whole asymmetric unit?

**Literature.** Carugo & Djinović-Carugo (J. Struct. Biol. 2012) counted molecules in actual contact with a reference molecule in PDB crystals. They found 3–18 contacts, about 10 on average, and fewer at higher symmetry. A Voronoi face is a geometric neighbour, not necessarily a contact.

**Method** (`pdb_chains.py`).
- Works for any space group and cell. The centroids of all protein chains (≥ 20 Cα) are expanded by the group operators (gemmi, via `research/p4/dvfast.Group`) and by lattice translations −2..2.
- The Voronoi tessellation of these points uses the real a, b, c, α, β, γ (scipy). Faces of each chain's cell = neighbouring chains.
- "ASU" mode puts one point at the centroid of all chains, as in `pdb_dv.py`.
- **Check:** for 35 single-chain P6₁22 entries the new code reproduces the old face counts exactly.
- **Sampling:** a random sample of 1 500 X-ray entries per group (seed 1). For I4₁32, all 158 entries.

**Result.**

| group | operations | entries | single-chain entries | ASU point: mean faces (max) | per chain: mean ± sd (max) | ASU point: odd face counts |
|---|---|---|---|---|---|---|
| P2₁ | 2 | 1 488 | 25% | 13.6 (16) | 15.3 ± 2.3 (29) | 0% |
| C2 | 4 | 1 476 | 33% | 13.5 (20) | 15.6 ± 2.6 (33) | 20% |
| P2₁2₁2₁ | 4 | 1 478 | 36% | 14.4 (20) | 15.3 ± 2.7 (27) | 0% |
| P3₂21 | 6 | 1 473 | 60% | 15.4 (23) | 15.9 ± 2.3 (26) | 41% |
| P6₁22 | 12 | 1 458 | 66% | 16.8 (24) | 16.8 ± 2.6 (29) | 50% |
| I4₁22 | 16 | 1 114 | 69% | 17.6 (29) | 17.1 ± 2.8 (31) | 57% |
| I4₁32 | 48 | 145 | 63% | 16.9 (20) | 16.6 ± 2.0 (23) | 71% |

The operations column counts the operations of the point group times the centring.

**Further measurements.**
- **Faces per chain grow with the number of chains in the ASU.** In P2₁: 13.6 for one chain, 14.7 for two, 15.2 for four, 15.9 for five or more. In P2₁2₁2₁: 14.5, 15.0, 15.6, 15.4.
- **Parity.** In P2₁ and P2₁2₁2₁ every single-point cell has an even number of faces. These groups contain no involution: no 2-fold rotation, mirror or inversion, since a 2₁ screw squared is a translation. So the neighbours x·g and x·g⁻¹ are always distinct, and faces come in pairs. Groups with 2-fold rotations (C2, P3₂21, P6₁22, I4₁22, I4₁32) have odd counts too.

**Conclusions.**
1. **Geometric neighbours increase with symmetry.** One point per molecule gives 13–14 faces in the low-symmetry groups (P2₁, C2, P2₁2₁2₁) and 17–18 in the high-symmetry ones. Carugo & Djinović-Carugo found the opposite trend for actual contacts. With more symmetry, more of the Voronoi neighbours are not touching the molecule.
2. **The difference between groups mostly disappears per chain.** With several chains in the ASU, the packing of chain centres is less regular, and every group gives 15–17 neighbours per chain. That is close to random packings (Poisson–Voronoi 15.5) and to the ~15 faces of a typical plant cell.
3. **Single molecules in P2₁ never exceed 16 faces** in the sample, and in P2₁2₁2₁ never exceed 20.

**Caveats.**
- **Chain centroids** stand in for whole chains. Long or elongated chains can contact neighbours that are not Voronoi neighbours of their centroids, and the other way round.
- **Sampling:** 1 500 random entries per group out of tens of thousands.
