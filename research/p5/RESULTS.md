# Problem 5: periodic tiling conjecture in dimension 3

**Question.** Must every finite set that tiles ℤ³ (or ℝ³) by translations also admit a fully periodic tiling?

**Status (18 Oct 2026).** Answered **no** by the OpenAI preprint "A translational tile with no fully periodic tiling in dimension three" (23 Sept 2026), together with a complete Lean 4 / Mathlib formalization. **We compiled and checked the formal proof ourselves:** it builds, uses only the standard axioms, and proves exactly the published challenge statement.

## Background (literature)
- **Dimension 1:** Newman (1977). Every tile of ℤ has a periodic tiling.
- **Dimension 2:** Bhattacharya (2020) proved the result for ℤ²; Greenfeld–Tao (2021) gave a quantitative version.
- **High dimension:** Greenfeld–Tao (Annals 2024, arXiv 2211.15847) disproved the conjecture.
- **Undecidability:** Greenfeld–Tao (JEMS 2025, arXiv 2309.09504) proved that translational monotiling of ℤᵈ, d ≥ 3, is algorithmically undecidable.
- **Correction.** Earlier versions of our site attributed this undecidability result to "Demaine–Langerman (Oct 2026)". That was wrong: no such paper exists. Fixed in `README.md` and `docs/js/problems.js`.

## The claim
- **Theorem 1.1 of the preprint.** A finite T ⊂ ℤ³ tiles ℤ³, but no tiling of ℤ³ by T is invariant under a subgroup of finite index.
- **Unit-cube thickening.** The same holds for T + [0,1]³ in ℝ³, with arbitrary real translations and tilings up to null sets.
- **Size of the tile.** The tile is astronomically large: the paper says no practical bound on the size is needed. It cannot be built or checked by computer search, so the only practical way to check the claim is the formal proof.

## What we checked

### Sources
- **OpenAI proof:** repository `openai/math`, commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`, directory `lean/OAI/Geometry/PeriodicTiling/`. It has 122 files, 16 899 lines and 886 theorems.
- **Toolchain:** Lean `v4.34.1`.
- **Mathlib:** commit `d13f23b723b8a846827a245b89c10fc7d3f11612`, exactly as in their `lake-manifest.json`.

### Checks
1. **Build.** `lake build OAI.Geometry.PeriodicTiling.Main` finished all 9049 jobs: the 122 proof modules plus Mathlib, which came from the official cache. There were no errors and no warnings in the proof files. Wall-clock time was 24 min on an M-series Mac (`lean_check/build.log`).
2. **Axioms.** `#print axioms OAI.PeriodicTilingThree.periodic_tiling_counterexample_and_minimality` gives `[propext, Classical.choice, Quot.sound]`. These are the three standard axioms of Lean/Mathlib. There is no `sorryAx` and no `Lean.ofReduceBool`, so no `native_decide`.
3. **Statement match.**
   - **Method.** `lean_check/LeanCheck/Dump.lean` dumps the theorem's type and every `OAI.*` constant it depends on, recursively, as raw kernel terms. It was run twice:
     - on the official challenge file (`ComparatorChallenges/PeriodicTilingThree.lean`, whose theorem is `sorry`);
     - on the solution environment.
   - **Coverage.** Both dumps contain the same 16 constants: `Lattice`, `Tiles`, `Period`, `FullyPeriodic`, `Space`, `castLattice`, `unitCube`, `Thickening`, `AETiles`, `latticeOfBasis`, `EuclideanFullyPeriodic` and `HasTileWithoutPeriodicComplement`, plus auxiliary proofs and the theorem itself.
   - **Result.** The dumps are identical except for automatically generated names of anonymous instance binders. Those names embed the source file name (`ComparatorChallenges.PeriodicTilingThree…` vs `OAI.Geometry.PeriodicTiling.Model…`). Binder names are irrelevant to the kernel, which uses de Bruijn indices.
   - **Mathlib.** Both environments use the same Mathlib commit, so the Mathlib constants in the statement are the same.
4. **No escape hatches.** The sources contain none of the following: `sorry`, `admit`, `axiom`, `opaque`, `unsafe`, `implemented_by`, `extern`, `native_decide`, `debug.*` options, or custom `elab`/`macro`/`syntax`/`#eval`/`run_cmd`.
5. **Kernel replay.** `leanchecker` (shipped with the toolchain) re-checks every declaration of each of the 122 compiled modules in the Lean kernel. Result: **122 of 122 modules pass, 0 failures** (`lean_check/leanchecker.log`).

### The formal statement, read by a human
The definitions are faithful to the mathematics.
- **`Tiles F A`:** (f, a) ↦ f + a is a bijection F × A → G, i.e. A ⊕ F = G.
- **`FullyPeriodic A`:** some finite-index subgroup consists of periods of A.
- **`Thickening T`:** the union of the unit cubes z + [0,1]ᵈ for z ∈ T.
- **`AETiles Ω A`:** almost every x lies in exactly one translate a + Ω, a ∈ A.
- **`EuclideanFullyPeriodic A`:** A is invariant under a full-rank lattice.

The theorem asserts two things:
- T exists and has all five properties: nonempty, tiles ℤ³, no periodic tiling of ℤ³, its thickening tiles ℝ³, and no tiling of ℝ³ by the thickening is periodic.
- 3 is the **least** dimension with a tile that has no periodic tiling. This part formalizes the d ≤ 2 results, including the planar periodicity theorem: see `LowerDimensions.lean` and `PlanarPeriodicity.lean`.

## Caveats
- **Scope.** A formal proof certifies the stated theorem, not the readability of the paper. The statement itself is the part a human has to check, and we did (above).
- **Trusted base.** We trust Lean's kernel, Mathlib and the compiler toolchain. That is standard.
- **What remains for people.** Independent human refereeing of the paper is still worth doing, but the result no longer depends on it.

## Reproduce
```
cd research/p5/lean_check
./fetch_sources.sh                      # OpenAI's 122 files into OAI/ (not committed)
lake exe cache get                      # Mathlib cache (~7.6 GB in .lake)
lake build OAI.Geometry.PeriodicTiling.Main ComparatorChallenges.PeriodicTilingThree LeanCheck.Dump
lake env lean CheckChallenge.lean       # writes dump_challenge.txt
lake env lean CheckSolution.lean        # writes dump_solution.txt, prints axioms
```
