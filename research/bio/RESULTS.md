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
