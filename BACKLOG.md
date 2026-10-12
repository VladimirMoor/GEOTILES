# Backlog

Ideas for later. For every item, the first step is to check the literature for existing solutions or work on the topic.

## Paused problems
- **Problem 3 (Heesch 7).** Verify Bašić's tile (`research/p3/tiles/basic6.txt`) with heesch-sat on a machine with 32–64 GB of RAM. Then look at families of Heesch-6 candidates.
- **Problem 5 (3D periodic tiling conjecture).** Independently verify the OpenAI counterexample (Sept 2026) before putting it in the catalog.

## Outside pure geometry (Voronoi cells and tilings in real data)
1. **Crystallography Open Database (COD).** Dirichlet cells of atoms and molecules in all space groups, compared with the record stereohedra of problem 4. *Done (11 Oct 2026): see `research/cod/RESULTS.md`.*
2. **3D cell tissues.** *Done (11 Oct 2026): see `research/tissue/RESULTS.md`; Weinberg types done (same symmetry and diversity as Voronoi of the centres at fixed F; excess of small cells). Still open: scutoids.* Use segmented 3D images of tissues (PlantSeg datasets, Arabidopsis ovule and meristem, Drosophila epithelia).
   - Measure the face-count distribution and the combinatorial types of cells.
   - Compare with Poisson–Voronoi (mean 15.54 faces), Kelvin (14) and proteins (~17).
   - Look for scutoids (Gómez-Gálvez et al. 2018).
3. **Virus capsids (VIPERdb).** *Checked 12 Oct 2026: thoroughly studied (Caspar–Klug, Mannige–Brooks, Twarock–Luque 2019). Covered as a review on the page `docs/nature.html`; a census of spherical Voronoi tilings of all PDB capsids would add little.* Tilings of the sphere: Caspar–Klug and Twarock's viral tiling theory. Classify the tilings that occur.
4. **Zeolites (IZA plus hypothetical databases).** *Checked 12 Oct 2026: done by others. Natural tilings of all frameworks: Anurova et al. 2010 (194), review 2018 (239), TTT tile collection. Isohedral simple tilings: Delgado-Friedrichs & O'Keeffe 2005. Covered as a short section on `docs/nature.html`.* Natural tilings of frameworks into cages (Blatov, ToposPro). Add tiles to the catalog and look for records.
5. **Foams and the Kelvin problem.** *Checked 12 Oct 2026: Weaire–Phelan is still the best; heavily studied. Side result: the most compact DV stereohedron in groups 75–230 is the truncated octahedron. Bezdek's bound, proved for parallelohedra in Sept 2026, appears to extend to stereohedra; see `research/p4/RESULTS.md`.* Is Weaire–Phelan optimal? Heavily studied (Kusner–Sullivan, Gabbrielli 2009). Low priority.
6. **Insect wings (dragonfly venation).** Test whether the wing cells are exactly Voronoi cells, using inverse-Voronoi recognition. 2D and already studied.
7. **Cosmic web (SDSS galaxies).** Voronoi face statistics against Poisson. Cosmology already uses Voronoi (ZOBOV, VIDE), so this is more a visual page than research.
8. **Columnar basalt and mud cracks.** 2D polygon statistics, well studied (Goehring). Low priority.

Not worth it as research: **honeycomb** and **wasp/hornet nests** (Hales 2001, Fejes Tóth 1964, Dauter–Jaskolski 2025; comb defects: Smith et al. 2023, Krishna et al. 2022). They are covered by the review page `docs/nature.html`, with our own bee-vs-Tóth computation (`research/nature/`). **Spider webs** are networks of threads, not tilings, and there is little data.

## Protein crystals (follow-ups)
- Split the asymmetric unit into chains instead of using its centroid.
- Other space groups: P2₁2₁2₁, P3₂21 and C2 are the most common.
