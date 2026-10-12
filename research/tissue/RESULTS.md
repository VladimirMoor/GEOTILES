# Cells in 3D plant tissues: face counts and Voronoi cells of their centres

**Question.** Cells in a tissue are convex-ish polyhedra that fill space. How many faces do they have? How close is the tissue to the Voronoi tessellation of its own cell centres, and to the Voronoi tessellation of random points (Poisson–Voronoi)?

**Literature (checked 11 Oct 2026).**
- 3D face and neighbour counts and the Euler, Lewis and Aboav–Weaire relations in a pseudostratified epithelium: Gómez, Iber et al., eLife 2021 (mouse lung).
- 3D connectivity networks of the shoot apical meristem: Jackson et al., Cell Systems 2019.
- Voronoi of nuclei as a model of epithelia: studied in **2D** (Sánchez-Gutiérrez et al., EMBO J 2015; fly and beetle embryos).
- Complete topology of 3D cells (Weinberg vectors): Lazar et al., PRL 2012, for grains and foams.
- We found no 3D test of "tissue vs Voronoi tessellation of the cell centres" on large segmented plant data.

**Data.** PlantSeg ground truth (Wolny et al., eLife 2020; OSF uzq3w).
- Arabidopsis ovules: 31 confocal volumes, voxel 0.235 × 0.150 × 0.150 µm (ds2x).
- Lateral root primordia: 26 light-sheet volumes, voxel 0.25 × 0.1625 × 0.1625 µm.
- Files were downloaded one at a time, processed and deleted.

**Method** (`cells.py`, `run_all.py`, `summary.py`).
- **Contacts.** Neighbouring voxels with different labels give a contact area for every pair of cells.
- **Faces.** A neighbour counts as a face if the contact area is at least a fraction *f* of the cell's surface. We report *f* = 0 (every contact) and *f* = 1%, 2%, 5% (noise removed).
- **Inner cells.** Only cells that touch neither background, nor the volume border, nor unlabelled regions.
- **Comparisons:**
  - the 3D Voronoi tessellation of the centroids of all cells of the same volume, with the same area threshold;
  - Poisson–Voronoi (random points) with the same threshold. Our Poisson–Voronoi gives 15.57 ± 3.33 faces, which matches the known values 15.54 and variance ≈ 11.

## Results (2026-10-11)

| | Ovules (17 827 inner cells) | Root primordia (325) | Voronoi of centres (ovules) | Poisson–Voronoi |
|---|---|---|---|---|
| faces, *f* = 0 | 15.16 ± 3.92 | 13.99 ± 3.85 | 14.71 ± 2.25 | 15.57 ± 3.33 |
| faces, *f* = 1% | 12.61 ± 2.84 | 11.63 ± 3.02 | 12.00 ± 1.86 | 11.95 ± 2.43 |
| same neighbours as Voronoi of centres (Jaccard) | 0.70 | 0.59 | – | – |

**Lewis law in 3D (ovules, *f* = 1%).** The mean volume of a cell with *n* faces, relative to the mean volume, rises steadily with *n*:

| n | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 |
|---|---|---|---|---|---|---|---|---|
| V/⟨V⟩ | 0.16 | 0.39 | 0.70 | 0.91 | 1.13 | 1.36 | 1.65 | 2.32 |

## Conclusions

1. **The mean face count of ovule cells (15.2) is close to Poisson–Voronoi (15.6) and to the Voronoi tessellation of the same centres (14.7).** Root primordia have fewer faces (14.0), near Kelvin's 14, but only 325 cells are complete.
2. **The distribution is much broader than for Voronoi of the same centres:** standard deviation 3.9 against 2.3, and even broader than random points (3.3). The tissue has many small cells with 4–9 faces and large cells with 20–30. The Voronoi tessellation of the centres evens out cell sizes; the tissue does not.
3. **Neighbour sets agree only partly.** On average 70% (ovules) or 59% (elongated root cells) of a cell's neighbours coincide with those of its Voronoi cell. In a 2D slice, 80% of the pixels fall in the same cell (see `docs/img/tissue-voronoi.png`). The tissue resembles Voronoi of its centres but is not one.
4. **Lewis law holds clearly in 3D:** volume grows monotonically with face count, from 0.16⟨V⟩ at 6 faces to 2.3⟨V⟩ at 20.

Next steps (see BACKLOG):
- classify the combinatorial types of cells with Weinberg vectors (Lazar et al.) and compare them with foams and grains;
- look for scutoids in curved epithelia.

## Complete topology of cells (Lazar et al. method), 2026-10-11

**Literature.** Lazar, Mason, MacPherson and Srolovitz (PRL 2012) classified grains and foam bubbles by Weinberg vectors. Grain growth gives symmetric cells "orders of magnitude" more often than Poisson–Voronoi. Lazar's research page lists biological tissue among the targets, but we found no published census of cell types in tissues.

**Method** (`topology.py`, `run_topo.py`, `topo_poisson.py`, `topo_summary.py`, `topo_byF.py`; check: `test_weinberg.py`).
- **Face graph from the voxels.** Faces are the neighbouring cells. Edges are triple lines: triples of labels in 2×2×2 blocks, seen in at least 2 blocks. Quadruple points add all four triples meeting at a vertex.
- **Slivers.** Faces with fewer than 3 edges are contracted onto an edge.
- **Accepted graphs.** Only 3-connected planar triangulations (the duals of simple polyhedra).
- **"Repaired" variant.** We additionally drop the weakest extra edges at nearly degenerate vertices.
- **Code.** For each accepted graph we compute the canonical Weinberg code, the symmetry order (reflections included) and the p-vector.
- **Validation.** The code gives the known symmetry orders for the tetrahedron (24), cube (48), triangular prism (12), truncated octahedron (48) and dodecahedron (120).
- **Comparisons:**
  - the Voronoi tessellation of the same centres, using the face graphs from the Delaunay triangulation;
  - Poisson–Voronoi, 21 149 inner cells.

**Coverage.**
- 31 ovules, 17 778 inner cells.
- Topology resolved for 8 989 cells (51%) in the strict variant and 10 963 (62%) with repair.
- Resolved cells have slightly fewer faces than all inner cells: 13.8 against 14.8. So all comparisons below are made **at a fixed number of faces F**.

**Overall numbers.** These raw shares mislead, because tissue has more small cells.

| | tissue | tissue (repaired) | Voronoi of same centres | Poisson–Voronoi |
|---|---|---|---|---|
| cells with non-trivial symmetry | 19.9% | 18.4% | 10.5% | 8.5% |
| same, reweighted to the tissue's F distribution | 19.9% | 19.4% | 19.2% | 16.6% |

**At fixed F** (share with non-trivial symmetry; number of distinct types among 150 random cells):

| F | tissue | Voronoi of same centres | Poisson–Voronoi |
|---|---|---|---|
| 10 | 0.58; 66 | 0.65; 47 | 0.50; 70 |
| 11 | 0.35; 91 | 0.39; 80 | 0.30; 111 |
| 12 | 0.24; 121 | 0.27; 112 | 0.19; 140 |
| 13 | 0.15; 142 | 0.15; 135 | 0.10; 149 |
| 14 | 0.08; 147 | 0.08; 145 | 0.06; 149 |
| 15 | 0.05; 149 | 0.04; 149 | 0.03; 149 |

**Most frequent types in tissue.** All of them are small cells, almost absent in both Voronoi models:

| F | p-vector (triangles, quads, pentagons, …) | symmetry | tissue | Voronoi of centres | Poisson |
|---|---|---|---|---|---|
| 6 | 2, 2, 2 | 4 | 1.1% | 0% | 0.03% |
| 5 | 2, 3 (triangular prism) | 12 | 0.8% | 0% | 0% |
| 7 | 1, 3, 3 | 6 | 0.6% | 0.02% | 0.08% |
| 4 | 4 (tetrahedron) | 24 | 0.3% | 0% | 0% |

**Conclusions.**
1. **At a fixed number of faces, tissue cells are as symmetric and as diverse as the Voronoi cells of the same centres.** Both are slightly more symmetric than Poisson–Voronoi. The raw excess of symmetric cells (20% against 10%) comes entirely from the many small cells in the tissue.
2. **Tissue does not behave like grain growth.** Lazar et al. found highly symmetric grains orders of magnitude more often than in Poisson–Voronoi; in tissue, at fixed F, the excess over Poisson is at most 1.5×, and it is the same for the Voronoi of the centres.
3. **The distinctive feature of tissue is a population of small cells with few faces:** tetrahedra, triangular prisms, 6- and 7-faced cells. By the Lewis-law data, cells with 4–6 faces have 7–16% of the mean volume. They are almost absent in both Voronoi models. They are present in the manually proofread ground truth, but we cannot exclude segmentation artefacts.
4. **The repaired and strict variants agree in every measure,** so the conclusions do not depend on the repair step.

## Tissue vs coarsening foam vs random cells: complete topology (2026-10-12)

**Question.** Is a plant tissue topologically like a foam, or like a random tessellation? Matzke compared tissue cells and foam bubbles by hand, using face counts and face types. Lazar et al. (PRL 2012) showed with Weinberg codes that coarsened grains are far more often symmetric than Poisson–Voronoi cells. Kraynik–Reinelt–van Swol (2003) compared simulated foam with Matzke's data. **We found no comparison of tissue and foam by complete cell topology.**

**Systems.** All four are processed by the same voxel pipeline (`topology.py`):
- tissue: 31 PlantSeg ovules;
- coarsening foam or grains: a 3D Potts model (`potts.py`), periodic box, 26-neighbour energy, kT = 0.5. It starts from a Voronoi tessellation of 40 000 seeds, and grains grow by curvature-driven coarsening, as in grain growth and foam ageing. Two runs:
  - 180³ box, snapshots at 400, 700 and 1000 steps;
  - 300³ box at tissue resolution, snapshots at 1400, 1800 and 2200 steps, i.e. 1027 → 489 grains of 26–55 thousand voxels;
- Poisson–Voronoi drawn in voxels at tissue resolution: 6 boxes of 300³ with 900 seeds each, about 30 000 voxels per cell;
- exact Poisson–Voronoi (Delaunay), for reference.

**Calibration of the pipeline** (voxelised Poisson–Voronoi against exact values: 15.5 faces, 8.5% symmetric):

| cell diameter (voxels) | resolved | mean F | symmetric |
|---|---|---|---|
| 12 | 31% | 13.4 | 17.3% |
| 17 | 48% | 14.1 | 15.4% |
| 24 | 60% | 14.6 | 12.4% |
| 31 (≈ tissue, which has a median of 31 000 voxels, ≈ 39 across) | 73% | 14.9 | 10.7% |

The pipeline loses small faces at coarse resolution, so cells look simpler and more symmetric. **This also biased our earlier tissue-vs-exact-Voronoi comparison.** Systems must be compared at equal resolution through the same pipeline, which is what the table below does.

**Result: share of cells with non-trivial symmetry at fixed F** (in brackets, distinct types among 100 random cells with that F):

| F | tissue | foam (Potts, 300³) | foam (Potts, 180³) | random, voxels | random, exact |
|---|---|---|---|---|---|
| 10 | 0.58 (47) | 0.82 (32) | 0.76 (25) | – | 0.50 (54) |
| 11 | 0.35 (76) | 0.55 (54) | 0.59 (35) | 0.20 (80) | 0.30 (75) |
| 12 | 0.24 (81) | 0.35 (74) | 0.41 (70) | 0.21 (97) | 0.19 (92) |
| 13 | 0.15 (94) | 0.18 (85) | 0.27 (79) | 0.11 (96) | 0.10 (96) |
| 14 | 0.08 (100) | – | 0.19 (99) | 0.06 (100) | 0.06 (100) |
| all cells | 0.199 | 0.448 | 0.470 | 0.115 | 0.085 |
| mean F ± sd | 13.8 ± 3.4 | 12.6 ± 4.9 | 12.5 ± 4.7 | 14.7 ± 3.2 | 15.5 ± 3.3 |

**Checks.**
- Potts at kT = 0.9 gives the same picture: 47% symmetric overall; at F = 11, 12 and 13: 0.60, 0.46 and 0.25. So it is not a temperature or roughening artefact.
- Perfect cubes, a possible lattice artefact, make up only 1–2% of the cells.

**Conclusions.**
1. **The coarsening foam (grain growth) is strongly ordered.** At fixed F it has 1.5–3 times more symmetric cells and fewer distinct types than random cells. Its face-count distribution is broad (sd 4.9). This agrees with Lazar et al.
2. **Plant tissue is much closer to the random tessellation than to the coarsening foam.**
   - At F ≥ 12 its symmetric share (0.24, 0.15, 0.08) lies between the voxel-random and foam values, but nearer to random.
   - Its type diversity is close to random.
   - Its face-count spread (3.4) is like random (3.2), not like foam (4.9).
3. Tissue is somewhat more ordered than random at small F, for example F = 11: 0.35 against 0.20. Cell division and wall mechanics leave a weaker imprint than surface-tension-driven coarsening.

**Caveats.**
- Only one tissue type (ovules).
- A Potts model stands in for foam; no real foam data was used. An open 40 GB raw foam tomography exists (Zenodo 15564270) but needs our own bubble segmentation.
- Topology is resolved for 51% of tissue cells against 73% of voxel-random cells at the same resolution. The tissue sample may be biased toward cleaner cells.
