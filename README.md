# GEOTILES

Open research on tilings. We work through a list of open problems about shapes that fill space. For each problem we first check the literature, then compute, preferably in exact arithmetic, and publish everything here.

**Site: https://geotiles-nine.vercel.app** (EN · RU · ES)

## Problems and results

| # | Problem | Status | Main results | Details |
|---|---|---|---|---|
| 1 | Equilateral convex polyhedra that tile space | in progress | **Theorem A**: every Minkowski sum of triangles and segments with total dimension ≤ 4 tiles. **Theorem B**: with total dimension ≥ 5 the tiling ones have measure zero. Two new equilateral space-fillers (9 and 10 faces). Classification of ≤ 8 faces, with 170 of 301 types rigorous. | [RESULTS](research/p1/RESULTS.md), [preprint](docs/paper/geotiles-minkowski-preprint.pdf) |
| 2 | New space-filling tetrahedra | in progress | The list 𝒜 of Chentouf–Sun (2023) has a misprint, and two of its "undecided" members are Sommerville tilers. None of the remaining 38 tiles face-to-face (12 of these exclusions are new). A golden box splits into 4 × #36 + #7. | [RESULTS](research/p2/RESULTS.md), [note](docs/paper/sporadic-tetrahedra-note.pdf) |
| 3 | A polygon with Heesch number 7 | paused | First Heesch census of polydrafters: 9.7 M shapes with ≤ 20 cells, maximum 3. Bašić's record tile reconstructed in exact coordinates (285 cells). Verifying it needs 32–64 GB of RAM. | [RESULTS](research/p3/RESULTS.md) |
| 4 | A space-filler with more than 38 faces | in progress | The bound 38 is not proved (Schmitt 2016 is a numerical search). Exact search over all groups 75–230 finds nothing above 38. New certified lower bounds for 10 groups. Engel's 38-hedron is provably rigid. | [RESULTS](research/p4/RESULTS.md), [note](docs/paper/dirichlet-stereohedra-note.pdf) |
| 5 | Periodic tiling conjecture in dimension 3 | resolved (preprint) | Counterexample by OpenAI (Sept 2026, not yet verified by people). Demaine–Langerman (Oct 2026) proved translational monotiling of ℝ³ undecidable. | — |
| 6–9 | Voronoi conjecture; full classification; decidability; boundedness of Heesch numbers | open | — | — |

Side project: Dirichlet cells of protein molecules in PDB crystals (`research/bio/`).

## The site

- **Catalog.** Equilateral space-fillers, record stereohedra (up to Engel's 38-hedron, with full space-group tilings) and, for comparison, solids that do not tile. Each comes with a 3D model, a tiling fragment, an explode slider, a slice view and automatic checks.
- **Lab.** A zonohedron builder with McMullen's criterion, and a polyhedron analyzer that runs the dihedral angle filter.
- **Research.** All problems with their status, results, preprints and notes.

## Structure

```
docs/                 static site (no build step)
  js/                 geometry, tiling, three.js viewer, i18n, catalog data
  paper/              compiled preprints and notes
research/p1/ … p4/    pipelines and RESULTS.md for problems 1–4
research/bio/         protein crystals from the PDB (site page: proteins.html)
research/paper*/      LaTeX sources of the preprint and notes
scripts/serve.py      local no-cache dev server
tools/                external tools (plantri, heesch-sat), git-ignored
```

## Run locally

```bash
python3 scripts/serve.py 8766
```

Then open http://localhost:8766.

## Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install numpy scipy mpmath sympy matplotlib gemmi
```

Each `research/pN/RESULTS.md` lists the scripts and commands for that problem. For example, the certified stereohedron bounds of problem 4 can be rechecked with:

```bash
cd research/p4 && ../../.venv/bin/python verify_certificates.py
```

## Deploy

`vercel.json` sets `docs/` as the output directory. Every push to `main` deploys.
