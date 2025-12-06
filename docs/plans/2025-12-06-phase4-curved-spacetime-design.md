# Phase 4: Curved Spacetime Sprinkling — Design Document

**Date:** 2025-12-06
**Status:** Approved, ready for implementation

---

## Overview

Extend the causal set simulator to sprinkle into curved spacetimes (starting with de Sitter) and verify that the Benincasa-Dowker action detects positive curvature.

### Goals

1. Create a `Spacetime<D>` trait abstraction for different geometries
2. Implement `Minkowski<D>` (refactor existing code)
3. Implement `DeSitter<D>` with conformal coordinates
4. Validate: S > 0 for de Sitter, S ≈ 0 for Minkowski

### Key Insight

For conformally flat spacetimes (Minkowski, de Sitter, FLRW), the causal structure is identical — only the volume measure differs. This makes the trait minimal:

```
Spacetime = Sampling + Causal Check + Curvature Info
```

---

## Architecture

### Core Trait

```rust
/// Trait defining a spacetime geometry for sprinkling
///
/// Coordinate convention:
/// - Index 0 is timelike
/// - Indices 1..D are spacelike
/// - Signature is (-,+,+,+)
pub trait Spacetime<const D: usize>: Clone + Send + Sync {
    /// Human-readable name for logging/testing
    fn name(&self) -> &'static str;

    /// Sample a point with correct volume weighting
    fn sample_point(&self, rng: &mut impl Rng) -> SVector<f64, D>;

    /// Check if p1 causally precedes p2
    fn causally_precedes(&self, p1: &SVector<f64, D>, p2: &SVector<f64, D>) -> bool;

    /// Ricci scalar curvature
    fn ricci_scalar(&self) -> f64;

    /// Proper volume of the sampling region
    fn volume(&self) -> f64;
}

/// Helper for conformally flat spacetimes
pub fn conformal_causal_check<const D: usize>(
    p1: &SVector<f64, D>,
    p2: &SVector<f64, D>
) -> bool {
    let dt = p2[0] - p1[0];
    if dt <= 0.0 { return false; }
    let dx_sq: f64 = (1..D).map(|i| (p2[i] - p1[i]).powi(2)).sum();
    dt * dt > dx_sq
}
```

### Minkowski Implementation

```rust
/// Flat Minkowski spacetime with sprinkling in a causal diamond
#[derive(Clone, Debug)]
pub struct Minkowski<const D: usize> {
    /// Half-height τ of the causal diamond (time spans [-τ, +τ])
    half_height: f64,
}

impl<const D: usize> Minkowski<D> {
    pub fn causal_diamond(half_height: f64) -> Self;
}

impl<const D: usize> Spacetime<D> for Minkowski<D> {
    fn name(&self) -> &'static str { "Minkowski-4D" /* etc */ }
    fn sample_point(&self, rng: &mut impl Rng) -> SVector<f64, D> { /* existing logic */ }
    fn causally_precedes(&self, p1, p2) -> bool { conformal_causal_check(p1, p2) }
    fn ricci_scalar(&self) -> f64 { 0.0 }
    fn volume(&self) -> f64 { /* C_d * τ^d */ }
}
```

### de Sitter Implementation

Conformal coordinates: `ds² = (1/H²η²)(-dη² + dx⃗²)` where η ∈ (-∞, 0)

```rust
/// de Sitter spacetime in conformal coordinates
#[derive(Clone, Debug)]
pub struct DeSitter<const D: usize> {
    h: f64,                    // Hubble parameter
    eta_range: (f64, f64),     // Conformal time bounds (both negative)
    spatial_half_width: f64,   // Spatial extent [-L, +L]
}

impl<const D: usize> Spacetime<D> for DeSitter<D> {
    fn name(&self) -> &'static str { "deSitter-4D" /* etc */ }

    fn sample_point(&self, rng: &mut impl Rng) -> SVector<f64, D> {
        // Importance sampling for η: √(-g) = 1/(H|η|)^D
        // Use inverse-CDF: η = -(a + u*(b-a))^(-1/(D-1))
        // Spatial: uniform (conformally flat)
    }

    fn causally_precedes(&self, p1, p2) -> bool {
        // Identical to Minkowski! Conformal invariance.
        conformal_causal_check(p1, p2)
    }

    fn ricci_scalar(&self) -> f64 {
        // R = d(d-1)H² → 12H² for 4D
        (D as f64) * (D as f64 - 1.0) * self.h * self.h
    }

    fn volume(&self) -> f64 {
        // ∫ √(-g) d^D x = (spatial vol) × ∫ dη/(H|η|)^D
    }
}
```

### Generic Sprinkler

```rust
pub struct Sprinkler<S, const D: usize>
where S: Spacetime<D>
{
    rng: ChaCha8Rng,
    spacetime: S,
    density: f64,
    // ...
}

impl<S, const D: usize> Sprinkler<S, D>
where S: Spacetime<D>
{
    pub fn new(spacetime: S, density: f64) -> Self;
    pub fn with_seed(spacetime: S, density: f64, seed: u64) -> Self;
    pub fn sprinkle(&mut self) -> SprinklingResult<D>;
    pub fn sprinkle_fixed(&mut self, n: usize) -> SprinklingResult<D>;
}
```

### CausalSet Construction

```rust
impl<const D: usize> CausalSet<D> {
    pub fn from_sprinkling<S: Spacetime<D>>(
        result: SprinklingResult<D>,
        spacetime: &S,
    ) -> Self {
        debug_assert_eq!(
            result.spacetime_name,
            spacetime.name(),
            "Spacetime mismatch"
        );
        // Build causal matrix using spacetime.causally_precedes()
    }
}
```

---

## Physics Reference

### de Sitter in Conformal Coordinates

- Metric: `ds² = (1/H²η²)(-dη² + dx⃗²)`
- Conformal time: η ∈ (-∞, 0), where η → 0⁻ is future infinity
- Volume element: `√(-g) = 1/(H|η|)^d`
- Ricci scalar: `R = d(d-1)H²` (constant, positive)
- Causal structure: **identical to Minkowski** (conformal invariance)

### Importance Sampling for η

For D dimensions, `√(-g) ∝ |η|^(-D)`:

```
CDF: P(η) ∝ ∫ dη/|η|^D = |η|^(1-D)/(1-D)

Inverse CDF:
  2D: η = -1/(1/|η_min| + u*(1/|η_max| - 1/|η_min|))
  3D: η = -(|η_min|^(-2) + u*(|η_max|^(-2) - |η_min|^(-2)))^(-1/2)
  4D: η = -(|η_min|^(-3) + u*(|η_max|^(-3) - |η_min|^(-3)))^(-1/3)

General: η = -(a + u*(b-a))^(-1/(D-1)) where a=|η_min|^(1-D), b=|η_max|^(1-D)
```

### Expected BD Action Behavior

| Spacetime | R | Expected S |
|-----------|---|------------|
| Minkowski | 0 | S ≈ 0 (fluctuations) |
| de Sitter | 12H² (4D) | S > 0 (positive curvature) |
| Anti-de Sitter | -12H² | S < 0 (negative curvature) |
| Schwarzschild | 0 (vacuum) | S ≈ 0 (Ricci-flat) |

---

## Validation Strategy

### Test 1: Sign Check (Most Critical)

```rust
#[test]
fn de_sitter_positive_action() {
    // Sprinkle multiple runs, compute t-statistic
    // Assert: mean_action / std_error > 2.0 (significantly positive)
}
```

### Test 2: Monotonicity

```rust
#[test]
fn de_sitter_action_scales_with_curvature() {
    // Test H = 0.05, 0.1, 0.2
    // Assert: S(H=0.2) > S(H=0.1) > S(H=0.05)
}
```

### Test 3: Flat vs Curved Comparison

```rust
#[test]
fn flat_vs_curved_action_comparison() {
    // Assert: |S_deSitter| >> |S_Minkowski|
}
```

### Test 4: Regression (Minkowski still works)

```rust
#[test]
fn minkowski_action_still_zero() {
    // Assert: action.is_consistent_with_flat()
}
```

### Statistical Robustness

Use t-statistic rather than simple mean > 0:
```rust
let t_statistic = mean_action / std_error;
assert!(t_statistic > 2.0, "Should be significantly positive");
```

---

## File Structure

```
crates/causet-core/src/
├── spacetime/           # NEW module
│   ├── mod.rs           # Spacetime trait + conformal_causal_check
│   ├── minkowski.rs     # Minkowski implementation
│   └── de_sitter.rs     # DeSitter implementation
├── sprinkling.rs        # MODIFIED: generic over Spacetime
├── causal_set.rs        # MODIFIED: from_sprinkling takes Spacetime
└── ...

tests/
└── phase4_curved_spacetime.rs  # NEW validation tests

examples/
└── de_sitter_validation.rs     # NEW validation example
```

---

## Migration Path

1. Create `spacetime/` module with trait
2. Implement `Minkowski<D>` (extract from existing code)
3. Make `Sprinkler` generic over `Spacetime`
4. Update `CausalSet::from_sprinkling` to take spacetime
5. Verify all existing tests pass
6. Implement `DeSitter<D>`
7. Add Phase 4 validation tests
8. Run validation, verify S > 0 for de Sitter

---

## Future Extensions

- **Anti-de Sitter**: Similar to de Sitter, expect S < 0
- **Schwarzschild**: Non-conformally-flat, needs geodesic integration for causal check
- **FLRW cosmologies**: Time-varying scale factor, conformally flat
- **Early termination optimization**: Skip causal checks for widely separated pairs

---

## References

- Benincasa-Dowker 2010 (arXiv:1001.2725) — BD action, de Sitter discussion
- Surya 2019 Living Reviews — Curved spacetime sprinkling
- Myrheim 1978 (CERN-TH-2538) — Ordering fraction values
