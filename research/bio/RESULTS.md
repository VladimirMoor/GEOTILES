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
