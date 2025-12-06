# Phase 3 Completion Report: Benincasa-Dowker Action

**Date:** 2025-01-05
**Status:** COMPLETE — BD action validated for flat spacetime

## Summary

Implemented the Benincasa-Dowker action, which is the discrete analog of the Einstein-Hilbert action. Verified that flat Minkowski spacetime produces S ≈ 0, indicating zero scalar curvature.

## Key Results

### Flat Spacetime Validation (N=500, 5 trials each)

| Dimension | Mean S | Std Dev | Result |
|-----------|--------|---------|--------|
| 2D (1+1)  | -61.6  | ±301.9  | ✓ Consistent with R=0 |
| 4D (3+1)  | -82.6  | ±360.8  | ✓ Consistent with R=0 |

Both dimensions show mean action close to zero with statistical fluctuations,
confirming flat spacetime has zero curvature as expected.

### Formula Verification (from BD 2010 primary source)

**2D action:**
```
S^(2) = N - 2N₁ + 4N₂ - 2N₃
```

**4D action:**
```
S^(4) = N - N₁ + 9N₂ - 16N₃ + 8N₄
```

Where:
- N = total elements
- N₁ = links (pairs with 0 elements strictly between)
- N₂ = pairs with 1 element strictly between
- N₃ = pairs with 2 elements strictly between
- N₄ = pairs with 3 elements strictly between

**Note:** BD paper uses (k+1)-element "inclusive" intervals, which equals our
k-element "strictly between" convention offset by 1.

## Implementation

### New File: `geometry/action.rs`

```rust
// 4D action (verified from BD 2010, arXiv:1001.2725v4)
pub fn bd_action_4d(counts: &IntervalCounts, n: usize) -> BDActionResult {
    let n1 = counts.n_k(0);  // links
    let n2 = counts.n_k(1);  // 1 between
    let n3 = counts.n_k(2);  // 2 between
    let n4 = counts.n_k(3);  // 3 between

    let action = n - n1 + 9*n2 - 16*n3 + 8*n4;
    // ...
}
```

### New CausalSet Methods

| Method | Description |
|--------|-------------|
| `bd_action_2d()` | Compute 2D BD action |
| `bd_action_4d()` | Compute 4D BD action |
| `bd_action()` | Dispatch based on dimension |

### New Types

- `BDActionResult` — Action value with metadata
- `BDCoefficients` — Dimension-specific coefficients

## Test Coverage

- **58 unit tests** (all passing)
- **7 doc tests** (all passing)
- **4 integration tests** (all passing)

New tests added:
- `test_bd_action_2d_formula` — Verify formula
- `test_bd_action_4d_formula` — Verify formula
- `test_bd_action_2d_flat_spacetime` — Validate S≈0
- `test_bd_action_4d_flat_spacetime` — Validate S≈0
- `test_bd_action_method_matches_dimension` — Dispatch test

## Physical Interpretation

The BD action encodes **scalar curvature** in discrete spacetime:

- **S ≈ 0**: Flat spacetime (Minkowski)
- **S > 0**: Positive curvature (e.g., de Sitter)
- **S < 0**: Negative curvature (e.g., Anti-de Sitter)

In the continuum limit:
```
lim S → ∫ √(-g) R d⁴x
```
This is the Einstein-Hilbert action that governs general relativity.

## Fluctuation Analysis

The fluctuations (std ≈ 300-360) are larger than the naive O(√N) ≈ 22 because:

1. The BD coefficients (±1, ±9, ±16, ±8) amplify statistical variations
2. N_k counts have correlated uncertainties
3. The true fluctuation scale is O(√N × coefficient_scale)

The important result is that **the mean is centered near zero**, confirming
flat spacetime has zero curvature.

## Reproducibility

```bash
# Run BD action validation
cargo run --example bd_action -p causet-core

# Run all tests
cargo test

# Run dimension validation (Phase 2)
cargo run --example dimension_test -p causet-core
```

## Next Steps: Phase 4

Curved spacetime sprinkling:
- de Sitter (expanding universe): expect S > 0
- Schwarzschild (black hole exterior): expect R ≠ 0 near horizon
- Compare action values to theoretical predictions

## References

- **Benincasa-Dowker (2010)** "Scalar Curvature of a Causal Set" arXiv:1001.2725v4
  - Line 599: 2D action formula
  - Lines 602-604: 4D action formula
  - Lines 605-606: N_k counting convention
