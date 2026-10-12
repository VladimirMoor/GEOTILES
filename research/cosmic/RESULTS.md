# Cosmic web: Voronoi cells of SDSS galaxies (2026-10-12)

**Literature.** Voronoi tessellations of galaxies are standard in cosmology:
- density estimates and void finders: ZOBOV, VIDE (review: arXiv:2012.08965);
- the volume distribution of cells and the number of neighbours against random points (Way et al., arXiv:1009.0387, SDSS DR7). Galaxy cells are far more varied in size.

We found no study of the complete combinatorial type of galaxy cells (Weinberg codes, symmetry at a fixed face count).

**Data** (`web.py`).
- SDSS DR17 main galaxy sample: r < 17.77 after extinction, spectroscopic redshift z = 0.02–0.10, sky region RA 130–235°, Dec −2–60°. 167k galaxies downloaded.
- Volume-limited cut M_r < −20.55 (complete to z = 0.1). This leaves 114 643 galaxies at a uniform mean density.
- Comoving coordinates for flat ΛCDM (H₀ = 70, Ω_m = 0.3).
- Exact 3D Voronoi cells (Delaunay), counted only in a core well inside the survey (RA 145–220°, Dec 8–50°, z 0.035–0.088): 38 209 cells.
- Baseline: the same number of uniform random points in the same region.

**Results.**

| | galaxies | random points |
|---|---|---|
| mean faces ± sd | 15.16 ± 5.04 | 15.53 ± 3.32 |
| spread of log cell volume | 1.39 | 0.45 |
| cells with non-trivial symmetry | 15.7% | 8.4% |

| faces F | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|
| symmetric, galaxies | 0.45 | 0.26 | 0.16 | 0.08 | 0.04 |
| symmetric, random | 0.51 | 0.29 | 0.19 | 0.09 | 0.05 |
| distinct types / 100, galaxies | 61 | 83 | 98 | 98 | 100 |
| distinct types / 100, random | 54 | 81 | 91 | 99 | 100 |

**Conclusions.**
1. Galaxy cells differ from random cells strongly in **size**: the volume spread is 3 times larger (filaments and voids). They also differ in the **spread of face counts** (sd 5.0 against 3.3).
2. But **at a fixed number of faces, their shapes are as random as those of random points**: the same symmetric share and the same variety of types. The higher overall symmetric share (15.7%) comes only from the many small, few-faced cells in dense regions.
3. So the "cosmic foam" is a foam in its geometry (sizes, walls, voids) but not in cell topology. A coarsening foam (Potts model, `research/tissue/RESULTS.md`) has 1.5–3 times more symmetric cells at the same F. In this sense galaxies behave like plant tissue: both are close to random.

**Caveats.**
- Redshift-space distortions ("fingers of God") stretch clusters along the line of sight.
- The survey footprint has holes that the random points do not reproduce.
- Both mainly affect cells near clusters and near edges.
- The core region and the volume-limited cut reduce, but do not remove, these effects.
