# Problem 4: a convex space-filler with more than 38 faces

## Literature check (2026-10-09)

**Record: 38 faces.** P. Engel, *Z. Kristallogr.* 154 (1981) 199–215, found a Dirichlet–Voronoi (DV) stereohedron, i.e. a plesiohedron, for the cubic group I4₁32 (IT 214). Its f-vector is (70, 106, 38). There are four combinatorial types with these numbers. 38 is the largest number of faces known for any convex monohedral tile in ℝ³.

**What is actually proved.** Wikipedia and MathWorld say that Schmitt (2016) proved 38 to be the maximum for plesiohedra. **That is not what the thesis says.**

M. Schmitt, *On Space Groups and Dirichlet–Voronoi Stereohedra*, PhD thesis, FU Berlin 2016, §2.3:
- It is a **numerical grid search** over the fundamental domains of the normalisers for 145 space groups: the tetragonal, trigonal, hexagonal and cubic ones. Triclinic, monoclinic and orthorhombic groups are omitted; they cannot produce more than 38 faces.
- It finds 3315 combinatorial types of DV-stereohedra. Every facet number from 4 to 38 occurs, and the maximum is attained only by Engel's examples.
- **It contains no proof that 38 is the maximum.**
- Schmitt himself calls problem (G) of Grünbaum–Shephard (1980) open: the least upper bound for the number of facets of a convex monohedral prototile in ℝⁿ. For n ≥ 3 it is not even known that the number is bounded.

**Rigorous upper bounds for DV-stereohedra**, proved group by group (Bochiş–Santos 2001, 2006; Sabariego–Santos 2008, 2011):

| Groups | Bound |
|---|---|
| Groups with reflections | ≤ 18 |
| Non-cubic groups | ≤ 80 |
| Full cubic groups | ≤ 25 |
| Quarter cubic groups | ≤ 92 (arXiv:0708.2114) |

So for plesiohedra the true maximum lies somewhere between 38 and 92.

**General bounds.**
- Every stereohedron (the tile of an isohedral tiling): Delone's bound 390, improved by Tarasov to 378.
- Non-isohedral monohedral tilings (anisohedral, k-isohedral): **no bound is known.**

**Not a record.** The 2026 "Josehedron" (arXiv:2604.07160) is a 12-face plesiohedron.

## Where room is left

1. **Plesiohedra, i.e. DV cells of one orbit.**
   - Schmitt's grid search covered this class thoroughly, numerically.
   - Possible gaps: very small parameter chambers that a grid misses, and special positions.
   - Low expected yield.
2. **Stereohedra that are not Dirichlet cells.**
   - These are convex tiles of isohedral tilings whose shape is not the Voronoi cell of an orbit. The gyrobifastigium is an example.
   - Bound 378; we know of **no systematic search**.
3. **Non-isohedral (k-isohedral) convex monohedral tilings.**
   - Example: Voronoi cells of several orbits that happen to be congruent.
   - **No bound and no systematic search known.**
