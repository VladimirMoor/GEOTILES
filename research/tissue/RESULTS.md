# Cells in 3D plant tissues: face counts and Voronoi cells of their centres

**Question.** Cells in a tissue are convex-ish polyhedra that fill space. How many faces do they have? How close is the tissue to the Voronoi tessellation of its own cell centres, and to the Voronoi tessellation of random points (Poisson–Voronoi)?

**Literature (checked 11 Oct 2026).**
- 3D face and neighbour counts and the Euler, Lewis and Aboav–Weaire relations in a pseudostratified epithelium: Gómez, Iber et al., eLife 2021 (mouse lung).
- 3D connectivity networks of the shoot apical meristem: Jackson et al., Cell Systems 2019.
- Voronoi of nuclei as a model of epithelia: studied in **2D** (Sánchez-Gutiérrez et al., EMBO J 2015; fly and beetle embryos).
- Complete topology of 3D cells (Weinberg vectors): Lazar et al., PRL 2012, for simulated grain growth (the model of coarsening grains and foams) and Poisson–Voronoi; no real foam.
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

**Literature.** Lazar, Mason, MacPherson and Srolovitz (PRL 2012) introduced Weinberg-vector classification for grains, bubbles and cells and applied it to simulated grain growth and to Poisson–Voronoi. Grain growth gives symmetric cells "orders of magnitude" more often than Poisson–Voronoi. Lazar's research page lists biological tissue among the targets, but we found no published census of cell types in tissues.

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
- A Potts model stands in for a coarsening foam. Real foam was added on 18 Oct 2026 (next section).
- Topology is resolved for 51% of tissue cells against 73% of voxel-random cells at the same resolution. The tissue sample may be biased toward cleaner cells.

## Real foam from X-ray tomography: complete bubble topology (2026-10-18)

**Question.** The comparison above used a simulated coarsening foam. What does the complete topology of real bubbles look like, and how far is tissue from it?

**Literature.**
- Matzke (1946) classified ~600 bubbles of a monodisperse soap foam by hand by face counts and face types (4-, 5-, 6-gons). The 13-faced 1-10-2 bubble was the most common.
- Real foams by optical tomography: Monnereau & Vignes-Adler (1998); Monnereau, Prunet-Foch & Vignes-Adler (2001). Mean face and edge counts are close to Matzke's; no complete topology.
- Simulated monodisperse foams: Kraynik, Reinelt & van Swol (PRE 2003, PRL 2004), mean F = 13.70. Sabik, Lutz & Evans (Soft Matter 2024) follow cell types during coarsening; in annealed monodisperse foam the 1-10-2 cell is ≈ 20% and the ten most common types ≈ 80%.
- Lazar et al. (PRL 2012) used no real foam, and **we found no Weinberg-type census of real foam bubbles from tomography.**

**Data.**
- Schott, Dollet, Santucci, Raufaste, Mokso, Nature Communications 2025 (tomo-rheoscopy), raw tomogram of series 2, step 102, Zenodo 15564270.
- Imaging: TOMCAT (Swiss Light Source), voxel 2.75 µm, 610 × 2016 × 2016 voxels.
- Foam: nearly monodisperse liquid foam (R ≈ 58 µm), liquid fraction ≈ 8%, between two rough plates (gap 1.5 mm) under very slow steady shear (shear rate ≤ 0.011 s⁻¹).
- We downloaded only slices z = 100–509 (13 GB, by HTTP range requests, `foam_fetch.py`) and stored them as uint8.

**Method.**
- `foam_seg.py`, as in FoamQuant:
  - Gaussian smoothing (σ = 1) and Otsu threshold for the liquid; liquid fraction found 9.6–10.2% per block.
  - Markers are maxima of the distance to the liquid (separation ≥ 12 voxels).
  - Watershed on −distance inside the sample disk (radius 940 voxels). The growth also fills the liquid, so the result is a space-filling tessellation whose triple lines are the Plateau borders.
- `foam_topo.py`: full resolution (bubble diameter ≈ 42 voxels, close to the tissue's ≈ 39), in 27 overlapping blocks. Each bubble is counted once, from the block whose core holds its centroid. Then the same pipeline as for tissue (`topology.py`). We also compute the exact Voronoi tessellation of all bubble centroids, for bubbles at least two diameters away from the ends of the range and from the sample edge.
- **Checks:**
  - Equivalent radius 59.1 µm, against 57.8 ± 1.0 µm in the paper.
  - Volume CV 0.18, including segmentation noise; the paper reports σ_R/R < 0.1.
  - Full resolution vs 2× binning on a 128 × 800 × 800 test block: mean F 13.80 vs 13.76, symmetric share 0.71 vs 0.76.
  - Visual check of label boundaries: `foam_test.py`, `docs/img/foam-real.png`.
- Commands: `foam_fetch.py 100 510`, `foam_topo.py 100 510`, `compare_systems.py`, `foam_summary.py` (→ `data/topo_foam/summary.json`), `foam_figure.py`.

**Result.** 23 047 inner bubbles; topology resolved for 18 093 (79%, against 51% for tissue).

| | real foam | Voronoi of its bubble centres | tissue | coarsening foam (Potts) | random (exact PV) |
|---|---|---|---|---|---|
| mean F ± sd | 13.88 ± 1.25 | 14.00 ± 1.31 | 13.84 ± 3.37 | 12.53 ± 4.69 | 15.53 ± 3.34 |
| symmetric, all cells | 0.68 | 0.56 | 0.20 | 0.47 | 0.085 |
| F = 11 | 0.71 (31) | 0.60 (32) | 0.35 (76) | 0.59 (46) | 0.30 (80) |
| F = 12 | **0.90 (13)** | 0.78 (38) | 0.24 (81) | 0.41 (66) | 0.19 (96) |
| F = 13 | **0.77 (21)** | 0.67 (43) | 0.15 (94) | 0.27 (79) | 0.10 (99) |
| F = 14 | **0.72 (40)** | 0.61 (48) | 0.08 (100) | 0.19 (94) | 0.06 (99) |
| F = 15 | 0.58 (50) | 0.47 (66) | 0.05 (100) | 0.11 (100) | 0.03 (100) |
| F = 16 | 0.30 (64) | 0.20 (77) | 0.04 (100) | 0.06 (99) | 0.01 (100) |

Each cell gives the share of cells with a non-trivial symmetry; in brackets, the number of distinct types among 100 random cells with that F.

**Real bubbles.**
- **Face counts:** F = 12: 8.7%, 13: 26.1%, 14: 34.3%, 15: 21.5%, 16: 6.4%. Cells with F ≤ 10 or F ≥ 18 make up 0.9%, mostly segmentation errors.
- **Most common types:**
  - 1-10-2 (the Matzke cell, symmetry order 4): 11.6%;
  - 1-10-3: 6.6%;
  - 2-8-4: 4.5%;
  - pentagonal dodecahedron 0-12-0 (order 120): 3.6%;
  - the 14-hedron 0-12-2 of the Weaire–Phelan foam (order 24): 3.5%;
  - 1-10-4: 3.5%.
- **Absent or nearly absent:**
  - Kelvin's truncated octahedron: 7 bubbles (0.04%);
  - cubes: 1 bubble.
- **Diversity:**
  - the ten most common types cover 43%;
  - 75% of bubbles have only 4-, 5- and 6-gonal faces.
- **Lewis law is flat:** the median volume of a bubble with F faces is 0.93 of the median at F = 11 and 1.03 at F = 17. In tissue it is 0.16 at F = 6 and 2.3 at F = 20.

**Bubbles vs the Voronoi cells of their own centres** (10 991 deep bubbles resolved in both):
- the same Weinberg type for 56% of bubbles, the same F for 73%;
- symmetric share 0.69 vs 0.60;
- at least one triangle 12.5% vs 21.9%, a heptagon or larger 21.1% vs 32.9%;
- the ten most common types cover 43% vs 31%.

**By distance from the axis** (shear rate grows linearly with radius): the symmetric share is 0.73, 0.71 and 0.65, and the ten most common types cover 50%, 47% and 40% in the inner, middle and outer thirds of the radius.

**Conclusions.**
1. **Real monodisperse foam is by far the most ordered system we have measured.**
   - At F = 12–14, 72–90% of bubbles have a symmetry, and only 13–40 distinct types occur among 100 bubbles.
   - That is 4–9 times more symmetric cells than in tissue, 5–12 times more than in random cells, and 2–4 times more than in the coarsening (polydisperse) Potts foam.
   - Matzke's 1-10-2 bubble is the most common type, as he found by hand in 1946.
2. **Most of this order is already in the packing of equal bubbles, and surface tension adds the rest.**
   - The Voronoi tessellation of the bubble centres is itself highly ordered (56% symmetric).
   - The real films remove a large part of the triangles and heptagons (Plateau's laws favour pentagons) and raise the symmetric share by 0.1 at every F.
   - A tissue, by contrast, is barely more ordered than random cells even though it is also a space-filling packing: its cells vary much more in size (Lewis law), and size disorder destroys topological order.
3. **Real foam is less ideal than simulated annealed foam:** the Matzke cell is 11.6% (≈ 20% in simulations), and the top ten types cover 43% (≈ 80%). Possible reasons are slow shear, segmentation noise and the 10% liquid fraction. Order decreases towards the rim, where the shear rate is highest, but image quality also drops there, so we do not attribute this to shear.
4. **For tissues, Matzke's analogy fails at the level of complete topology.** Plant cells are not "bubbles": topologically they are much closer to random Voronoi cells than to a real foam.

**Caveats.**
- **One foam:** one sample, one time step; a wet (≈ 10%), sheared, nearly monodisperse foam.
- **Liquid assigned to bubbles:** the tessellation splits the liquid between bubbles, and faces are bubble adjacencies in that tessellation. Near-contacts across a thick Plateau border can add or remove small faces.
- **Unresolved bubbles:** 21% of bubbles fail the 3-connected planarity test, mostly near-degenerate vertices; the repaired variant is stored as `cells_repaired`.
- **Resolution bias:** at this resolution voxelised cells look slightly more symmetric than they are (calibration above). The bias is of the order of +0.02 to +0.05 and does not affect the comparison, since the bubbles and tissue cells span the same number of voxels.

## A second real foam: open-cell polyurethane (2026-10-18)

**Question.** The liquid foam above is nearly monodisperse, and most of its order is already in the packing of equal bubbles. Does a polydisperse real foam also show topological order, and where does it come from?

**Literature.**
- Montminy, Tannenbaum & Macosko (J. Colloid Interface Sci. 2004) reconstructed six polyurethane foams from tomography: ~8 000 struts, 4 000 windows, but only **376 cells**. They report window shapes (55% pentagons against ~70% in aqueous foams) and cell face counts, not complete cell types.
- Jang, Kraynik & Kyriakides (IJSS 2008): micro-CT of open-cell PU and aluminium foams; cells are irregular polyhedra with 9–17 faces.
- Kraynik, Reinelt & van Swol (PRL 2004), simulations: mean F decreases with polydispersity, and quadrilaterals take over from pentagons.
- We found no Weinberg-type census of a solid foam either.

**Data.**
- DaRUS doi:10.18419/darus-3010 (CC BY 4.0, Univ. Stuttgart): micro-XRCT of an open-cell reticulated PUR foam, 10 pores per inch, cylinder 52.4 × 44.3 mm, undeformed state t0.
- We used the authors' own segmented stack (strut skeleton against void, 593 slices of 751 × 751, voxel ≈ 0.07 mm, 6 MB).

**Method.**
- `pur_seg.py`: the cells have no membranes, only struts along the edges, so cells are separated where the void narrows at the windows. Distance to the struts (Gaussian σ = 1), markers at h-maxima, watershed inside the cylinder (radius 360 voxels). Struts are absorbed into cells, so they become triple lines.
- **Choice of h** (`pur_check.py`, central block of 300 × 400 × 400):

| h | triple-line voxels within 2 voxels of a strut | strut voxels within 3 voxels of a triple line | cells |
|---|---|---|---|
| 0.5 | 0.948 | 0.859 | 706 |
| 1 | 0.966 | 0.843 | 630 |
| 2 | 0.976 | 0.811 | 556 |
| 4 | 0.985 | 0.737 | 444 |
| 8 | 0.984 | 0.483 | 192 |

  - In a real open-cell foam every cell edge is a strut. So the first column measures false edges (over-segmentation) and the second measures missed edges (merged cells).
  - We use h = 1; h = 2 serves as a check.
- `pur_topo.py`: 27 overlapping blocks, as for the liquid foam. Cell equivalent diameter ≈ 58 voxels (4.4 mm), which is finer than the tissue's resolution.
- **Cut surface excluded.** Cells within ~1.2 cell diameters of the sample's cut surface (r > 290 voxels) have abnormally few faces: 9.7 against 13.5 inside, at normal volume. Struts there are cut or damaged, so we exclude them. Along the axis there is no such trend.

**Result.** 1 641 inner cells, 1 396 resolved; 1 103 of them are away from the cut surface.

| | PU foam | Voronoi of its cell centres | real liquid foam | coarsening foam (Potts) | tissue | random (exact PV) |
|---|---|---|---|---|---|---|
| mean F ± sd | 13.51 ± 3.03 | 14.83 ± 2.08 | 13.88 ± 1.25 | 12.53 ± 4.69 | 13.84 ± 3.37 | 15.53 ± 3.34 |
| symmetric, all cells | 0.46 | 0.13 | 0.68 | 0.47 | 0.20 | 0.085 |
| F = 12 | 0.61 (44) | – | 0.90 (13) | 0.41 (72) | 0.24 (81) | 0.19 (96) |
| F = 13 | 0.46 (47) | – | 0.77 (21) | 0.27 (78) | 0.15 (94) | 0.10 (99) |
| F = 14 | 0.38 (65) | 0.13 (96) | 0.72 (40) | 0.19 (94) | 0.08 (100) | 0.06 (99) |
| F = 15 | 0.29 (78) | 0.07 (98) | 0.58 (50) | 0.11 (97) | 0.05 (100) | 0.03 (100) |
| F = 16 | 0.18 (86) | 0.06 (99) | 0.30 (64) | 0.06 (100) | 0.04 (100) | 0.01 (100) |

As before, each cell gives the share of cells with a non-trivial symmetry, and in brackets the number of distinct types among 100 random cells with that F. A dash means fewer than 100 cells.

**Further measurements.**
- **Check with h = 2:** 1 077 cells, mean F 13.53 ± 2.87. The symmetric share is higher by 0.04–0.07 at every F (F = 14: 0.44), so the comparison is unchanged.
- **Polydispersity:**
  - volume CV 0.48;
  - Lewis law: the median volume at F = 11 is 0.61 of the overall median, at F = 17 it is 1.42 (liquid foam: 0.93–1.03; tissue: 0.16 at F = 6 and 2.3 at F = 20).
- **Face shapes:** 27% of cells have a triangular face and 48% have a face with 7 or more sides (liquid foam: 12.5% and 21%). The most common types are 4-4-0 (symmetry order 8), 3-6-0 (order 12), 2-8-1 and 3-6-1, each 2–3%. The Matzke 1-10-2 cell is only 1.5–2%.

**Conclusions.**
1. **A polydisperse solid foam is also strongly ordered.**
   - At F = 12–16 it has 1.5–3 times more symmetric cells than the coarsening Potts foam and 2.5–6 times more than tissue.
   - It has far fewer distinct types than tissue or random cells.
2. **Here the order comes from the films, not from the packing.**
   - The Voronoi cells of the PU cell centres are almost as disordered as random cells (F = 14: 0.13 against 0.06).
   - For the monodisperse liquid foam the Voronoi of centres was already ordered (0.61). Together the two foams separate the two sources of order.
   - Equal sizes give an ordered packing. Surface tension (Plateau's laws, which the PU foam inherited from its liquid state) gives order on top of any packing.
3. **Tissue has neither.** Its face-count spread is about the same as the PU foam's (sd 3.4 against 3.0), and its cells vary even more in size, but its symmetric shares are near random. Tissue walls are not shaped like soap films in topology.
4. **Face statistics follow Kraynik et al. (2004).** Polydispersity lowers the share of pentagons and of the Matzke cell and raises quadrilateral-rich types (4-4-0, 3-6-x).

**Caveats.**
- **Small sample:** 1 103 cells, one sample. Shares at fixed F have a statistical error of about ±0.04.
- **Segmentation:** the open-cell segmentation misses ~16% of strut length as edges, merging some faces.
- **Rise anisotropy:** PU cells are elongated along the rise direction. This does not change topology, but the Voronoi of centroids is computed without correcting for it.
- **Different system:** the foam is a solid replica of a liquid foam, reticulated (membranes removed), so it is a different system from the liquid foam.

## Real metal grains: pure iron from diffraction contrast tomography (2026-10-24)

**Question.** Lazar et al. (PRL 2012) found that simulated grain growth (mean curvature flow) produces highly symmetric grains far more often than random cells, and our Potts model agrees. Do real grains show this?

**Literature.**
- Complete-topology censuses of grains exist only for simulations: Lazar et al. 2012, and later grain-growth simulation studies. We found none for a measured 3D polycrystal.
- Experiments: Bhattacharya et al. (Science 2021) showed with HEDM that real grain-boundary velocities in Ni are not correlated with curvature, i.e. real grain growth is not simple curvature flow.

**Data.**
- Ånes, Heinig, Wang, Oddershede, … Zhang, Lauridsen, Zenodo 20323210 (CC BY 4.0): laboratory DCT of pure α-iron during grain growth (sample AGG100), state T0, GrainMapper3D HDF5.
- 3292 × 371 × 370 voxels of 3 µm, a rod of 1.1 × 1.1 × 9.9 mm.
- The grain map fills the sample: 60 135 grains, median equivalent diameter 15.9 voxels (48 µm), volume CV 1.9.
- The archive was fetched in 16 parallel range requests, because a single stream ran at ≈ 35 kB/s.

**Method.** `fe_topo.py`.
- GrainId is used directly as the label volume, with 0 outside the sample mask, and the same pipeline as for tissue (`topology.py`). Grains touching the mask or the ends of the volume are excluded.
- Exact Voronoi of all grain centroids, for grains deeper than 32 voxels.
- **Resolution matters** (calibration above): grains here are ≈ 16 voxels across. We therefore compare with systems drawn at the same resolution:
  - the Potts snapshot at 400 steps (180³ box, ≈ 2 000 grains, ≈ 18 voxels);
  - voxelised Poisson–Voronoi with cell diameter ≈ 14 voxels (`pv_N180_n4000`).

**Result.** 41 647 inner grains; topology resolved for 23 037 (55%; voxelised random cells at the same resolution: 31%).

| F | iron grains | Voronoi of iron grain centres | grain growth, Potts t400 (same res.) | Potts t400–1000 | random, voxels D ≈ 14 | tissue | random (exact) |
|---|---|---|---|---|---|---|---|
| mean F ± sd | 11.38 ± 4.64 | 14.83 ± 3.25 | 12.58 ± 4.78 | 12.53 ± 4.69 | 13.44 ± 2.90 | 13.84 ± 3.37 | 15.53 ± 3.34 |
| 10 | 0.53 | 0.71 | 0.75 | 0.76 | 0.62 | 0.58 | 0.50 |
| 11 | **0.31** | 0.47 | 0.54 | 0.59 | 0.29 | 0.35 | 0.30 |
| 12 | **0.21** | 0.32 | 0.35 | 0.41 | 0.18 | 0.24 | 0.19 |
| 13 | **0.12** | 0.16 | 0.28 | 0.27 | 0.04 | 0.15 | 0.10 |
| 14 | **0.06** | 0.10 | 0.18 | 0.19 | 0.04 | 0.08 | 0.06 |
| 15 | **0.03** | 0.07 | 0.17 | 0.11 | 0.02 | 0.05 | 0.03 |

Each cell gives the share of cells with a non-trivial symmetry. Distinct types among 100 iron grains: 75, 88, 98 and 99 at F = 11–14, against 43, 69, 85 and 97 in the Potts model.

**Further measurements.**
- **Face counts are grain-growth-like:** mean 11.4, sd 4.6. That is broad, like the Potts model (12.5 ± 4.7) and unlike random cells (15.5 ± 3.3). 17% of grains have ≤ 7 faces, and 1.7% are tetrahedra.
- **Lewis law:** median volume relative to the overall median is 0.82 at F = 11 and 2.34 at F = 17.
- **Most common types:** small cells such as the 6-faced 2-2-2 (4.5%), the triangular prism (3.4%), 1-3-3-1 and 1-3-3 (2.5% each).
- **Bubble types are absent:** grains with only 4-, 5- and 6-gonal faces are 11%, against 75% in the liquid foam.
- **Robustness to reconstruction quality** (`fe_check.py`, z = 1000–1600, 3 936 resolved grains):
  - grains in the top and bottom halves of mean LabDCT completeness have the same symmetric shares at F = 12: 0.24 and 0.25 (F = 11: 0.37 and 0.30);
  - large and small grains also agree;
  - a stricter triple-line threshold cannot be used: at 16 voxels per grain real edges are only a few voxels long, and MIN_TRIPLE = 6 drops mean F to 6.

**Conclusions.**
1. **Real iron grains are topologically as disordered as random cells at fixed F.** Their symmetric shares lie on the voxel-random and tissue curves and are 1.7–6 times below simulated grain growth at the same resolution. They are even below the Voronoi cells of their own centres.
2. **Face counts follow grain growth, but grain shapes do not.** The size-related statistics (broad F distribution, Lewis law, many small few-faced grains) look like grain growth. The topological order that mean curvature flow produces (Lazar et al.) is absent.
3. **This is consistent with Bhattacharya et al. (2021):** real boundaries do not move by curvature alone. Anisotropic boundary energy and mobility, solute drag and pinning plausibly break the order that isotropic curvature flow creates.
4. **Overall order of systems at fixed F:** monodisperse liquid foam > polydisperse PU foam > simulated grain growth (Potts) > real iron grains ≈ plant tissue ≈ random cells.
   - Surface-tension-dominated cells are ordered.
   - Real solids and living tissue are not.

**Caveats.**
- **One sample, one state:** T0 of a series prepared for abnormal grain growth.
- **Resolution:** coarse, ≈ 16 voxels per grain, though compared at matched resolution.
- **Reconstruction:** LabDCT reconstruction may miss grains below a few voxels and smooths boundaries. A bias common to all grains cannot be excluded from one dataset.
- **Next steps:** later states T1–T4 of the same sample, and a higher-resolution grain map (HEDM, 3D EBSD).

## Scutoids in the plant epidermis (2026-10-25)

**Question.** In curved animal epithelia many cells are scutoids: their neighbours on the apical surface differ from those on the basal surface (Gómez-Gálvez et al., Nat Commun 2018). Plant cells are glued by shared walls and cannot slide past each other. Are there scutoids in a plant epidermis, and how many?

**Literature.**
- Gómez-Gálvez et al. 2018 report the share of scutoids by tissue:

| tissue | scutoids | surface ratio |
|---|---|---|
| Drosophila salivary gland | 75 ± 17% | 7 ± 3 |
| Drosophila embryo, folds | 50 ± 15% | 1.6 |
| Drosophila embryo, outside folds | 15 ± 4% | 1.2 |
| egg chamber, stage 4 / stage 8 | 20 ± 8% / 10 ± 3% | |
| zebrafish epiboly | 2.9 ± 1.5% | |

- In their Voronoi tube model the share grows with the surface ratio.
- **We found no study of scutoids in plants.**

**Data.** The same 31 PlantSeg ovule volumes, re-downloaded from OSF one at a time and deleted after processing. Four parallel downloads were needed, since OSF ran at about 200 kB/s per connection.

**Method.** `scutoid.py`, `scutoid_all.py`, `scutoid_summary.py`, `scutoid_figure.py`.
- **The epithelium.** Cells touching the background that also touch cells which never touch it, i.e. cells spanning the outer layer.
  - A first version took every cell touching the background and found ~40% "scutoids". Slices showed this was an artefact: in ovules, cells are often stacked two deep near the surface, and the integuments are sheets with background on both sides.
- **Relative depth** in the layer: t = d_out / (d_out + d_in), with distances in µm to the outside and to the tissue below the layer.
- **Wall bands.** For each pair of layer cells, the wall area in the outer band (t < 0.2) and in the inner band (t > 0.8).
- **Scutoid criterion.** A cell is a scutoid when it takes part in an apico-basal T1 exchange:
  - a pair (a, b) has a wall ≥ A µm² in one band and no wall at all in the other;
  - and a pair (c, d) of common neighbours of a and b does the opposite.
  - A one-sided change without the compensating pair is reported separately; it is mostly a wall that thins out near a surface.
- **Cells judged:** spanning cells with a basal surface, not touching the volume border or ignore regions. 3 317 cells in 31 ovules; median layer thickness 5.6 µm.
- **Surface ratio** (outer area / inner area of the judged cells): median 1.16 (quartiles 1.01–1.27).
- **Visual check:** a verified exchange in N_405_A, `docs/img/scutoid-ovule.png`. Cells A and B touch on the outer surface, C and D on the inner.
- **Two null models on the same layer and the same cells** (`column_model`): each cell becomes a straight column, and at each depth level (20 levels) a layer voxel goes to the nearest column axis, i.e. a Voronoi tessellation at each depth, like the Voronoi tubes of Gómez-Gálvez et al.
  - "normal": axis along the outward surface normal at the centre of the cell's outer face;
  - "real axes": axis through the real centres of the cell's outer and inner faces, which keeps each cell's real tilt.

**Result** (share of judged cells in at least one apico-basal T1 exchange):

| wall threshold A | real epidermis | Voronoi columns, real axes | Voronoi columns, normals | ovules with real < model (real axes) |
|---|---|---|---|---|
| 0.5 µm² | 12.8% | 33.3% | 31.3% | 28 of 31 |
| 1.0 µm² | **9.0%** | 23.5% | 21.4% | 26 of 31 |
| 2.0 µm² | 4.5% | 12.5% | 11.4% | 26 of 31 |

**Conclusions.**
1. **Scutoids exist in plants.** A clean apico-basal neighbour exchange is seen in the ovule epidermis, and about 9% of epidermal cells (4.5–13% depending on the threshold) take part in one.
2. **They are 2.5–2.8 times rarer than geometry alone predicts.** The same cells, as Voronoi columns with the same tilt, would exchange neighbours 2.5–2.8 times more often. This holds in 26–28 of 31 ovules.
   - Plant cells cannot intercalate. Their walls are laid down by division and then only stretch, so it is plausible that they keep the same neighbours through the layer more than a free-boundary (Voronoi) packing would.
3. **Comparison with animals.** At a similar surface ratio (1.16 against 1.2) the ovule epidermis has 9% scutoids against 15 ± 4% in the Drosophila embryo outside folds. The methods differ, so this is only indicative.

**Caveats.**
- **Few judged cells:** only spanning epidermal cells away from the borders and ignore regions count: 3 317 of 14 378 spanning cells. Some ovules have only a few dozen judged cells.
- **Segmentation errors** can create or destroy T1 exchanges. The models have none, and errors more likely add exchanges to the real layer than remove them, so if anything the deficit of real scutoids is underestimated.
- **Voronoi columns are a reference, not a mechanical model.** Walls in a Voronoi column model are planar, and part of the difference may come from the curved walls of real cells.
