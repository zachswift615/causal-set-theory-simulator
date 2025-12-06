# Session Handoff: Phases 2 & 3 Complete

**Date:** 2025-01-05
**Session Focus:** Geometric observables and Benincasa-Dowker action
**Status:** Phases 2 & 3 complete, ready for Phase 4

---

## What We Accomplished

### Phase 2: Myrheim-Meyer Dimension Estimator

**Goal:** Extract spacetime dimension from pure causal structure.

**Key Discovery:** The gamma function formula often cited in literature is WRONG:
```
f(d) = Γ(d+1) · Γ(d/2) / [4 · Γ(3d/2)]  ← Does NOT reproduce Myrheim's values!
```

**Solution:** Use Myrheim 1978 Table I values directly with log-linear interpolation:

| Spacetime | Ordering Fraction f |
|-----------|---------------------|
| 2D (1+1)  | 1/2 = 0.500         |
| 3D (2+1)  | 8/35 ≈ 0.229        |
| 4D (3+1)  | 1/10 = 0.100        |

**Validation Results (N=500):**
- 2D: d = 2.01 ± 0.01 (0.5% error)
- 3D: d = 3.04 ± 0.09 (1.3% error)
- 4D: d = 4.07 ± 0.03 (1.7% error)

**Implementation:** `geometry/dimension.rs`
- Log-linear interpolation through anchor points
- `DimensionEstimate` struct with extrapolation flag for d > 4
- `DimensionStatistics` with uncertainty δd ≈ d/√N

---

### Phase 3: Benincasa-Dowker Action

**Goal:** Implement the discrete Einstein-Hilbert action.

**Formulas (verified from BD 2010, arXiv:1001.2725v4):**
```
2D: S^(2) = N - 2N₁ + 4N₂ - 2N₃
4D: S^(4) = N - N₁ + 9N₂ - 16N₃ + 8N₄
```

**N_k Counting Convention:**
- N₁ = links (pairs with 0 elements strictly between) = our `counts[0]`
- N₂ = pairs with 1 element strictly between = our `counts[1]`
- N₃ = pairs with 2 elements strictly between = our `counts[2]`
- N₄ = pairs with 3 elements strictly between = our `counts[3]`

**Note:** BD paper uses (k+1)-element "inclusive" intervals, offset by 1 from our convention.

**Validation Results (N=500, flat Minkowski):**
- 2D: S = -61.6 ± 301.9 (consistent with R=0)
- 4D: S = -82.6 ± 360.8 (consistent with R=0)

**Fluctuation Context from BD Paper:**
- N=5000: μ = 9.35, σ = 134.8
- N=20000: μ = 1.12, σ = 58.8
- Our results are consistent with this scaling

**Implementation:** `geometry/action.rs`
- `bd_action_2d()`, `bd_action_4d()` functions
- `BDActionResult` with normalized action S/√N
- `is_consistent_with_flat()` helper

---

## Code Structure After This Session

```
causet-core/src/
├── lib.rs                 # Updated with geometry exports
├── point.rs               # SpacetimePoint<D>
├── sprinkling.rs          # Poisson process, CausalDiamond
├── relations.rs           # CausalMatrix, LinkMatrix
├── causal_set.rs          # CausalSet main struct
└── geometry/              # NEW in this session
    ├── mod.rs             # Module root + CausalSet extensions
    ├── dimension.rs       # Myrheim-Meyer estimator
    ├── chains.rs          # Longest chain / proper time
    ├── intervals.rs       # N_k counting
    └── action.rs          # Benincasa-Dowker action
```

**New CausalSet Methods:**
- `estimated_dimension()` → `DimensionEstimate`
- `dimension_statistics()` → `DimensionStatistics`
- `longest_chain(x, y)` → `ChainResult`
- `interval_counts()` → `IntervalCounts`
- `bd_action_2d()` → `BDActionResult`
- `bd_action_4d()` → `BDActionResult`
- `bd_action()` → dispatches based on D

---

## Key Lessons Learned

### 1. Always Verify Against Primary Sources

Both the Myrheim-Meyer formula AND the BD coefficients needed verification:
- Myrheim-Meyer: Secondary sources had wrong values (0.424, 0.333)
- BD action: Our α/β parameterization had inconsistency; direct formula is correct

### 2. Sprinkling Method Matters

`sprinkle_fixed()` gives correct sampling distribution.
`sprinkle_diamond()` with density can produce biased results for small diamonds.

Always use:
```rust
let diamond = CausalDiamond::<D>::symmetric(10.0);
let mut sprinkler = Sprinkler::with_seed(seed, 1.0);
let result = sprinkler.sprinkle_fixed(&diamond, n_points);
```

### 3. Counting Convention

Our `interval_counts().n_k(k)` = pairs with k elements **strictly between**.
BD paper's N_k = (k+1)-element **inclusive** intervals.
So: BD's N₁ = our `n_k(0)`, BD's N₂ = our `n_k(1)`, etc.

---

## Git History

```
541b4f0 Phase 1: Core sprinkling engine with validated Myrheim-Meyer ordering fractions
0943e7c Phase 2: Geometric observables (Myrheim-Meyer dimension estimator)
0a1d612 Phase 3: Benincasa-Dowker action (discrete Einstein-Hilbert)
```

---

## Test Coverage

- 58 unit tests
- 7 doc tests
- 4 integration tests

All passing.

---

## Phase 4: Curved Spacetime Sprinkling

### Goal
Sprinkle into curved spacetimes and verify non-zero BD action.

### Test Cases (from claude.ai analysis)

| Spacetime | Expected S | Reason |
|-----------|-----------|--------|
| de Sitter | S > 0 | Positive constant curvature R = 12H² |
| Anti-de Sitter | S < 0 | Negative constant curvature |
| Schwarzschild | S ≈ 0 | Ricci-flat (R = 0 in vacuum!) |

**Important insight:** Schwarzschild is Ricci-flat, so BD action should still be ~0.
The curvature there is in the Weyl tensor, not Ricci scalar.

**De Sitter is the best test case** — constant positive R throughout.

### Implementation Needs

1. **de Sitter sprinkling** with proper volume element
2. **Coordinate systems** that work for sprinkling (global de Sitter, EF for Schwarzschild)
3. **Causal relation tests** in curved spacetime
4. **Volume element** corrections for non-flat metrics

### References for Phase 4

- BD 2010 (arXiv:1001.2725v4) — has de Sitter discussion
- Surya 2019 Living Reviews — curved spacetime sprinkling
- Glaser papers — generalized coefficients for any dimension

---

## Commands to Verify Current State

```bash
# Build
cargo build

# Run all tests
cargo test

# Phase 2 validation (dimension recovery)
cargo run --example dimension_test -p causet-core

# Phase 3 validation (BD action)
cargo run --example bd_action -p causet-core

# Phase 1 validation (ordering fractions)
cargo run --example basic_sprinkling -p causet-core
```

---

## Files Changed This Session

```
Created:
  crates/causet-core/src/geometry/mod.rs
  crates/causet-core/src/geometry/dimension.rs
  crates/causet-core/src/geometry/chains.rs
  crates/causet-core/src/geometry/intervals.rs
  crates/causet-core/src/geometry/action.rs
  crates/causet-core/examples/dimension_test.rs
  crates/causet-core/examples/bd_action.rs
  docs/plans/2025-01-05-phase2-geometric-observables-design.md
  PHASE2_COMPLETION.md
  PHASE3_COMPLETION.md
  docs/SESSION_HANDOFF_PHASES_2_3.md (this file)

Modified:
  crates/causet-core/src/lib.rs
  CAUSET_PROJECT_PROMPT.md
```

---

## Summary

Phases 2 and 3 establish that:

1. **Spacetime dimension emerges from causal structure** (Myrheim-Meyer)
2. **Scalar curvature is encoded in interval counting** (Benincasa-Dowker)
3. **Flat spacetime correctly gives S ≈ 0**

Phase 4 will test curved spacetimes to see if we can detect curvature through the BD action.
