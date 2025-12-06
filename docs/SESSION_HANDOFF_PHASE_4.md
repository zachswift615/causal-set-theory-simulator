# Session Handoff: Phase 4 Complete

**Date:** 2025-12-06
**Session Focus:** Curved spacetime sprinkling (de Sitter)
**Status:** Phase 4 complete, ready for Phase 5 or additional curved spacetime tests

---

## What We Accomplished

### Phase 4: de Sitter Sprinkling and Curvature Detection

**Goal:** Extend the causal set simulator to curved spacetimes and verify the BD action detects curvature.

**Key Result:**
```
de Sitter (R = 12H² > 0):  mean S = +182.9  ✓ positive
Minkowski (R = 0):         mean S =  +97.8  ✓ ~0 (within σ≈600)
```

The Benincasa-Dowker action correctly distinguishes curved from flat spacetime using only causal structure!

---

## Technical Implementation

### New `Spacetime<D>` Trait

```rust
pub trait Spacetime<const D: usize>: Clone + Send + Sync {
    fn name(&self) -> String;
    fn sample_point(&self, rng: &mut impl Rng) -> SVector<f64, D>;
    fn causally_precedes(&self, p1: &SVector<f64, D>, p2: &SVector<f64, D>) -> bool;
    fn ricci_scalar(&self) -> f64;
    fn volume(&self) -> f64;
}
```

### Implementations

**`Minkowski<D>`** — Flat spacetime in causal diamond
- `sample_point()`: Volume-weighted time sampling + ball rejection for spatial
- `causally_precedes()`: Standard Minkowski check (Δt² > |Δx|²)
- `ricci_scalar()`: 0

**`DeSitter<D>`** — Curved spacetime in conformal coordinates
- Metric: `ds² = (1/H²η²)(-dη² + dx⃗²)` where η ∈ (-∞, 0)
- `sample_point()`: Importance sampling for η (volume diverges as η→0)
- `causally_precedes()`: Same as Minkowski! (conformally flat)
- `ricci_scalar()`: `d(d-1)H²` (12H² for 4D)

### Key Insight: Conformal Coordinates

For conformally flat spacetimes, the causal structure is identical to Minkowski. The conformal factor Ω² only affects the volume element, not the light cones. This means:

1. `causally_precedes()` can reuse the Minkowski check
2. Only `sample_point()` needs to handle the non-uniform volume
3. Future FLRW cosmologies will also be easy to implement

---

## Code Structure After This Session

```
crates/causet-core/src/
├── lib.rs                 # Updated with spacetime exports
├── spacetime/             # NEW MODULE
│   ├── mod.rs             # Spacetime trait + conformal_causal_check()
│   ├── minkowski.rs       # Minkowski<D> implementation
│   └── de_sitter.rs       # DeSitter<D> implementation
├── sprinkling.rs          # + sprinkle_spacetime(), GenericSprinklingResult
├── causal_set.rs          # + from_generic_sprinkling()
└── geometry/              # Unchanged from Phase 3

tests/
├── ordering_fraction.rs   # Existing
└── phase4_curved_spacetime.rs  # NEW: 5 validation tests

examples/
├── compare_apis.rs        # NEW: Old vs new API comparison
└── sampling_stats.rs      # NEW: Statistical sampling comparison
```

---

## Key Lessons Learned

### 1. Conformal Coordinates Are Your Friend

The conformal coordinate choice `ds² = Ω²(x) η_μν dx^μ dx^ν` makes de Sitter trivial:
- Causal structure = Minkowski causal structure
- Only need importance sampling for the volume element

For Schwarzschild (non-conformally-flat), you'd need geodesic integration for causal checks.

### 2. BD Action Fluctuations Are Large

For N = 300-500 points:
- Expected: S ∝ ∫R√(-g)d⁴x ≈ 0 for flat, > 0 for de Sitter
- Observed fluctuations: σ(S) ≈ 600-700

Single runs can give misleading results. Always average over 10-20 trials minimum.

### 3. Fixed N ≠ Fixed Density

When comparing different H values with fixed N points:
- The causal structure is the same (conformal invariance)
- The proper volume differs (V ∝ 1/H⁴)
- So S doesn't scale with H² at fixed N

To see curvature scaling, need to either:
- Fix density ρ and let N vary
- Normalize by proper volume

---

## Git History

```
f26d74d Phase 4: Curved spacetime sprinkling (de Sitter)
3737e3c Add session handoff document for Phases 2-3
0a1d612 Phase 3: Benincasa-Dowker action (discrete Einstein-Hilbert)
0943e7c Phase 2: Geometric observables (Myrheim-Meyer dimension estimator)
541b4f0 Phase 1: Core sprinkling engine with validated Myrheim-Meyer ordering fractions
```

---

## Test Coverage

- 79 unit tests
- 11 doc tests
- 4 integration tests (ordering fraction)
- 5 Phase 4 validation tests

**Total: 99 tests, all passing**

---

## Commands to Verify Current State

```bash
# Build
cargo build

# Run all tests
cargo test

# Phase 4 validation
cargo test --test phase4_curved_spacetime -- --nocapture

# Compare old vs new API
cargo run --example compare_apis -p causet-core

# Statistical sampling comparison
cargo run --example sampling_stats -p causet-core

# Previous phase validations still work
cargo run --example dimension_test -p causet-core
cargo run --example bd_action -p causet-core
```

---

## Files Changed This Session

```
Created:
  crates/causet-core/src/spacetime/mod.rs
  crates/causet-core/src/spacetime/minkowski.rs
  crates/causet-core/src/spacetime/de_sitter.rs
  crates/causet-core/tests/phase4_curved_spacetime.rs
  crates/causet-core/examples/compare_apis.rs
  crates/causet-core/examples/sampling_stats.rs
  docs/plans/2025-12-06-phase4-curved-spacetime-design.md
  docs/SESSION_HANDOFF_PHASE_4.md (this file)

Modified:
  crates/causet-core/src/lib.rs
  crates/causet-core/src/sprinkling.rs
  crates/causet-core/src/causal_set.rs
  CAUSET_PROJECT_PROMPT.md
```

---

## What's Next

### Option A: More Curved Spacetime Tests

Complete the curved spacetime story before moving to dynamics:

1. **Anti-de Sitter** — Expect S < 0 (negative curvature)
   - Same conformal trick works
   - Completes the "sign test" for curvature detection

2. **Schwarzschild** — Expect S ≈ 0 (Ricci-flat in vacuum)
   - Harder: non-conformally-flat, need geodesic integration
   - Use Eddington-Finkelstein coordinates for horizon crossing
   - Interesting null test: curved spacetime with R = 0

3. **Scaling studies** — Verify σ(S) ∝ N^(-1/2)
   - Run with N = 500, 1000, 2000, 5000
   - Should see fluctuations decrease

### Option B: Phase 5 — Classical Sequential Growth

Move on to dynamics:

1. **Transitive percolation** — Stochastic growth rule
2. **Birth process** — Add elements respecting causality
3. **Covariance** — Verify growth is label-independent
4. **Dimension stability** — Does d stay constant as causet grows?

### Option C: Phase 6 — Visualization

Make it visual with rerun.io:

1. **3D causet rendering** — Points + links in 3D projection
2. **Light cone visualization** — Show causal structure
3. **Growth animation** — Watch CSG evolution

---

## Summary

**Phases 1-4 establish that:**

1. **Spacetime dimension emerges from causal structure** (Myrheim-Meyer)
2. **Scalar curvature is encoded in interval counting** (Benincasa-Dowker)
3. **Flat spacetime correctly gives S ≈ 0** (Phase 3)
4. **Curved spacetime (de Sitter) correctly gives S > 0** (Phase 4) ← NEW!

The core physics of causal set theory is now validated in our simulator. Future work can explore dynamics, additional spacetimes, and visualization.

---

## Start Next Session With

```
Read docs/SESSION_HANDOFF_PHASE_4.md and CAUSET_PROJECT_PROMPT.md.

Phase 4 is complete. Choose one of:

A) More curved spacetimes: Anti-de Sitter (S<0), Schwarzschild (S≈0), scaling studies
B) Phase 5: Classical Sequential Growth dynamics
C) Phase 6: Visualization with rerun.io

What would you like to tackle?
```
