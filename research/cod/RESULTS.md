# COD: Dirichlet cells of atom orbits in real crystals

**Question.** Take one atom in a crystal and the orbit of its position under the space group. The Voronoi cell of that orbit is a convex polyhedron whose copies fill space: a Dirichlet stereohedron of the same kind as the problem-4 records (Engel's 38). How many faces do these cells have in real crystals? Do real structures come near the record of their group?

**Literature (checked 10 Oct 2026).**
- Blatov, Serezhkin and coauthors have large VDP statistics over ICSD and CSD; the 14-face cell is the most frequent. Those cells are built from **all** atoms of the structure, or from one element's sublattice. They are not single-orbit stereohedra, and they are not compared with the record bounds of each group.
- We found no study of single-orbit cells against per-group records.

**Method** (`cod_dv.py`, `summary.py`).
- Take all COD entries in 14 non-centrosymmetric groups where the record stereohedra are large.
- Read the cell, symmetry operations and atom sites from the CIF with gemmi.
- Skip an entry if its CIF operations differ from the standard setting (9 entries) or if its metric does not fit.
- Exclude H and D atoms.
- For every remaining site, build the DV cell of its orbit with `research/p4/dvfast.py`.
- **Baseline.** Uniformly random points x, at c/a values drawn from the same COD entries (3000 points per group).
- The top 3 general-position sites of every group were re-certified in exact rational arithmetic (`research/p4/exact.py`). All 42 agree with the float counts.

## Results (2026-10-11)

"Max" is the maximum over general-position sites with full occupancy; the bracket gives the maximum over all sites, disordered ones included. "Exp." is the number of general sites with ≥ max faces expected for random points: count × the random frequency. A dash means the frequency is below 1/3000 and was not resolved.

| Group | Record | Entries | General sites | Max (all) | ≥ max obs. / exp. | Mean faces COD / random |
|---|---|---|---|---|---|---|
| I4₁32 | 38 | 46 | 374 | 28 (28) | 2 / 1.3 | 17.25 / 17.19 |
| I4₁22 | 35 | 133 | 2425 | 29 (30) | 1 / 0.8 | 17.39 / 17.51 |
| P6₁22 | 34 | 157 | 4370 | 30 (30) | 3 / 5.8 | 16.69 / 16.63 |
| P6₅22 | 34 | 145 | 4078 | 31 (31) | 3 / 4.1 | 16.82 / 16.72 |
| I-42d | 33 | 501 | 9101 | 28 (28) | 2 / – | 16.76 / 16.92 |
| P4₃32 | 29 | 92 | 1275 | 25 (26) | 1 / 1.3 | 15.39 / 15.69 |
| P4₁32 | 29 | 91 | 804 | 25 (25) | 1 / 2.1 | 15.32 / 15.62 |
| I4₁cd | 24 | 188 | 5679 | 22 (22) | 1 / 3.8 | 15.48 / 15.50 |
| I2₁3 | 24 | 107 | 1905 | **24** (24) | 1 / 1.9 | 14.39 / 14.84 |
| P2₁3 | 24 | 596 | 10283 | **24** (24) | 6 / 13.7 | 12.79 / 12.86 |
| P4₁2₁2 | 29 | 964 | 33477 | 26 (26) | 4 / – | 16.16 / 16.13 |
| P4₃2₁2 | 29 | 852 | 30083 | 26 (26) | 4 / – | 16.07 / 16.06 |
| I4₁ | 26 | 135 | 7085 | **26** (26) | 38 / 71 | 18.08 / 17.96 |
| P3₁21 | 25 | 546 | 14104 | 24 (24) | 6 / 4.7 | 15.55 / 15.66 |

The record for each group is Schmitt's 2016 value, or our certified value for I-42d (33) and I4₁cd (24).

## Conclusions

1. **Atom positions behave like random points for this statistic.** The mean face counts agree with the random baseline to within 0.3 in every group, and the tail counts (≥ max) agree in order of magnitude. There is no sign that crystals prefer or avoid many-faced orbit cells.
2. **The group record is reached in three groups**, I2₁3, P2₁3 and I4₁, with exact certificates. Examples are COD 7223002 (C10, now in the site catalog), 7707836 and 7243061. In each case the random frequency of the record is large enough that such hits are expected. In I4₁ about 1% of random points reach 26.
3. In groups with high records (I4₁32: 38, I4₁22: 35, P6₁22: 34) the record cells are very rare for random points, and crystals stay at 74–91% of the record.
4. **Parity.** In P2₁3 every general position has an even face count. Neighbours g·x and g⁻¹·x pair up, and the pairing has no fixed point because the group has no involutions: its order-2 elements are 2₁ screws, whose square is a translation. This is elementary and surely known; it explains the "every other bin" histograms.

Data: COD (Gražulis et al., Nucleic Acids Res. 40 (2012) D420–D427). Per-entry results are in `data/cod_*.json`, the summary in `data/summary.json`. The CIF cache (1.2 GB) is git-ignored.
