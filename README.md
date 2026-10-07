# GEOTILES

Open research on tilings, starting with one question: **which convex polyhedra with all edges equal fill space?**

**Site: https://geotiles-nine.vercel.app**

The site has three languages (EN, RU, ES) and these parts:

- **Catalog.** Equilateral space-fillers and, for comparison, solids that do not tile. Each comes with a 3D model, a tiling fragment, an explode slider, a slice view and automatic checks.
- **Lab.** A zonohedron builder that applies McMullen's criterion and constructs the translation tiling. An analyzer that takes arbitrary vertices and runs the dihedral angle filter.
- **Research.** Literature review, plan and log.

## Verification built into the site

- **Dihedral angle filter.** This is a necessary condition for a monotiling: every dihedral angle must be part of a sum Σ kᵢαᵢ = 360° or 180°.
- **Volume check.** V(tile) / V(fundamental domain) must equal 1.
- **Coverage check.** 1500 random points in a ball covered by the fragment must each lie in exactly one tile.

## Structure

```
docs/            static site (deployed as is, no build step)
  js/geom.js     convex hull, dihedral angles, zonohedra, McMullen criterion, angle filter
  js/tiling.js   tiling fragments and coverage check
  js/viewer.js   three.js viewer
  js/i18n.js     translations
  js/catalog-data.js  catalog entries
  vendor/three/  three.js r170 (MIT)
scripts/serve.py local no-cache dev server
problem1_literature.md  literature review for problem 1 (RU)
```

## Run locally

```bash
python3 scripts/serve.py 8766
```

Then open http://localhost:8766.

## Deploy

The site is static. `vercel.json` sets `docs/` as the output directory, so a Vercel project imported from this repository needs no further settings.
