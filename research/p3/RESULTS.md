# Problem 3: a polygon with Heesch number 7

Started 2026-10-09.

## Literature check (2026-10-09)

**Record.** The largest known finite Heesch number is **6**. No tile with Heesch number 7 or more has been published.

| Year | Result | Source |
|---|---|---|
| 2021 | First tile with Heesch number 6. It is a 75-edge polygon on the kisrhombille grid, i.e. a polydrafter. | B. Bašić, "A Figure with Heesch Number 6", *Math. Intelligencer* 43 (2021) 50–53, doi:10.1007/s00283-020-10034-w |
| 2019–2026 | Further tiles with Heesch number 6, found by Dave Smith and confirmed by Bašić and Bakker. They are reported as strict H = 6, with a 7th corona refuted. | arXiv:2610.02579; https://heesch-atlas.netlify.app |
| 2001 / 2004 | Earlier record of 5, held by Mann's family. Mann and Thomas reached 5 with edge-marked polyforms. | https://faculty.washington.edu/cemann/Heesch.pdf |

**Kaplan's exhaustive census of unmarked polyforms** (*Contrib. Discrete Math.* 17(2) 2022; https://cs.uwaterloo.ca/~csk/heesch/). Hc is the Heesch number when the outermost corona must have no holes; Hh allows holes.

| Polyform | Sizes covered | Maximum Hc |
|---|---|---|
| Polyominoes | n ≤ 19 | 3 |
| Polyhexes | n ≤ 17 | 4 |
| Polyiamonds | n ≤ 24 | 4 |
| Polykites | about n ≤ 16 | about 2–3 |

**No published Heesch census of polydrafters**, even though every known H = 6 tile appears to live on the drafter (kisrhombille) grid.

**Tool: heesch-sat** (https://github.com/isohedral/heesch-sat, C. S. Kaplan, BSD-3). It is a SAT-based computation of Heesch numbers for many polyform grids, including `-drafter`.

## Setup

**Build** (macOS, Homebrew). heesch-sat lives in `tools/heesch-sat`, which is git-ignored.
- Dependencies: CryptoMiniSat 5.16 and boost from Homebrew.
- Compile with `-std=c++20`. The lambdas in the templates need C++20.
- Link `isohedral.o` into `sat`. The upstream Makefile omits it.

**Validation.** Running on the free 6-hexes reproduces the README exactly:
- 81 shapes;
- 1 inconclusive (the 2-anisohedral one);
- 4 non-tilers: three with Hc = 1 and one with Hc = 2;
- 76 isohedral tilers.

## Plan

1. **Census of free polydrafters by size** (`census.sh`).
   - The fast pass stops at 5 coronas.
   - A tile with H ≥ 7 would show up as "inconclusive" at that level, so every inconclusive shape and every shape with Hc = 5 is re-run individually with a higher `-maxlevel` and a time limit.
   - Shape counts grow by a factor of about 2.4 per cell: 5129 free 12-drafters.
2. Report the maximum Heesch number of polydrafters as a function of size. This is new data.
3. Later: parametric families around the known H = 6 tiles (Bašić, Smith), which have about 140 drafters each, far beyond an exhaustive census.

## Census of free polydrafters, sizes 2–20 (2026-10-09)

**Method** (`census.sh`, `recheck.py`).
- **Fast pass.** heesch-sat with `-isohedral -maxlevel 5`.
- **Re-check.** Every inconclusive shape and every shape with Hc ≥ 3 is re-run individually with `-periodic -maxlevel 8/9`. The `-periodic` switch is undocumented; it searches for periodic, non-isohedral tilings.
- Reports: `data/drafterN_report.txt`. Re-check logs: `data/recheck_*.log`. The shape files themselves are large and git-ignored.

**Total:** about 9.7 million free polydrafters with up to 20 cells, about 45 minutes of CPU on a laptop. In the table, a dash means no shapes of that kind.

| n | shapes | Hc=0 | Hc=1 | Hc=2 | Hc=3 | isohedral | inconclusive → resolved |
|---|---|---|---|---|---|---|---|
| 2 | 3 | – | – | – | – | 3 | – |
| 3 | 3 | 1 | – | – | – | 2 | – |
| 4 | 9 | 4 | – | – | – | 5 | – |
| 5 | 14 | 10 | 3 | 1 | – | – | – |
| 6 | 38 | 11 | – | – | – | 27 | – |
| 7 | 73 | 64 | 5 | 2 | – | – | 2 → 2-anisohedral |
| 8 | 178 | 138 | 8 | – | – | 29 | 3 → 2-anisohedral |
| 9 | 388 | 368 | 4 | 1 | – | 10 | 5 → 2-anisohedral |
| 10 | 933 | 888 | 37 | 1 | – | – | 7 → 2-anisohedral |
| 11 | 2 139 | 2 119 | 18 | 2 | – | – | – |
| 12 | 5 129 | 4 281 | 24 | – | – | 816 | 8 → 2-anisohedral |
| 13 | 12 096 | 12 057 | 34 | 1 | **1** | – | 3 → 2-anisohedral |
| 14 | 29 105 | 28 407 | 127 | 3 | **1** | 561 | 6 → 2-anisohedral |
| 15 | 69 594 | 69 418 | 135 | 3 | – | 16 | 22 → 2-anisohedral |
| 16 | 167 715 | 166 772 | 482 | 7 | – | 413 | 41 → 2-anisohedral (one needs a low maxlevel) |
| 17 | 403 251 | 403 086 | 162 | 2 | **1** | – | – |
| 18 | 972 551 | 959 169 | 516 | – | – | 12 866 | – |
| 19 | 2 343 836 | 2 343 536 | 294 | 3 | **1** | – | 2 → 2-anisohedral |
| 20 | 5 658 732 | 5 656 446 | 2 183 | 21 | **1** | – | 71 → 2-anisohedral; 10 see below |

**Result.** The maximum hole-free Heesch number of a polydrafter with at most 20 cells is **Hc = 3**. It is attained by exactly one shape at each of n = 13, 14, 17, 19, 20, and all five values were confirmed in the re-check.

**The 10 unresolved 20-drafters.**
- Their `-periodic` check crashes with a segfault inside `periodic.h`, after "Hmm 1".
- The crash happens only after 8 coronas have been built. So each of these shapes is either a tiler or has H ≥ 8. In particular, none can have Heesch number exactly 7.
- All 10 are long, thin strips with parallel long sides, and they visibly tile by stacking.

**Interpretation.**
- Small polydrafters behave like the other polyform families, whose maxima are 3–4.
- The known Heesch-6 tiles are much larger: Bašić's tile has about 140 drafters.
- Exhaustive search will not reach that size; the census grows by a factor of about 2.4 per cell.
- The next step has to be targeted families built from the known H = 6 tiles.

## Bašić's Heesch-6 tile in heesch-sat coordinates (2026-10-10)

**Source.** Wikimedia Commons, `File:A_polygon_with_Heesch_number_6.svg`, a CorelDRAW drawing of Bašić's tile with its 6 coronas.
- The drawing has 338 tile copies in total: the central tile and the coronas, with 6, 10, 13, 32, 52 and 55 tiles.
- Each copy is drawn as 6 regular hexagons and 69 small triangles of shape 30-60-90, i.e. drafters.

**In drafter units.**
- One equilateral triangle of the kisrhombille is 6 drafters, and one hexagon is 36 drafters.
- So the tile has 6·36 + 69 = **285 drafters**. This matches the description "6 hexes + 11.5 triangles", since 69 = 11.5 · 6.

**Conversion** (`svg_to_drafter.py`).
- Each hexagon is subdivided into drafters.
- A similarity is fitted to heesch-sat's drafter grid. Its rotation comes from all edge directions, and its scale and translation from a scale scan plus ICP. The fitted SVG triangle side is 30.57 units against 7 grid units.
- One copy (class `fil1`) aligns with maximum centroid error 0.13 cell, and all 285 cells are distinct.
- The central copy does not align: it has 3 malformed triangles in the SVG, one with angles 26.6°/90°/63.4°.
- Result: `tiles/basic6.txt`, in heesch-sat format.

**Running heesch-sat** (`guarded_sat.py`, with a memory watchdog polled every 0.5 s).
- Even at `-maxlevel 2`, with or without `-isohedral`, memory passes **4.4 GB within 4 s**, before any corona is solved.
- The cost is in the set-up. All placements of a neighbouring copy are enumerated: 12 orientations × translations × 285 cells, about 10⁶ placements.
- A first unguarded attempt was killed by the OS at more than 2 GB after 85 s.

**Conclusion.** Reproducing H = 6 for Bašić's tile with heesch-sat needs a machine with substantially more RAM, probably 32–64 GB, or a memory-lean reimplementation of the neighbour-placement step. That rules it out on this 16 GB laptop. Bakker reports similar resource limits (CryptoMiniSat clause ceilings) for tiles of this size.
