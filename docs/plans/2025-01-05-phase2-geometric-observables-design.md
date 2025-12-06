# Phase 2 Design: Geometric Observables

**Date:** 2025-01-05
**Status:** Approved
**Milestone:** Sprinkle into 2D, 3D, 4D diamonds → recover correct dimensions within error bars

## Overview

Phase 2 implements geometric observables that extract spacetime geometry from pure causal structure — the core claim of causal set theory. The key deliverable is the Myrheim-Meyer dimension estimator, validated against sprinkled causets.

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Crate organization | Extend `causet-core` | Fast iteration > premature modularity |
| Interpolation method | Log-linear | 4 anchor points don't justify PCHIP complexity |
| d > 4 handling | Struct with extrapolation flag | Research needs to know when extrapolating |

## Module Structure

```
causet-core/src/
├── lib.rs              # Add geometry module export
├── point.rs            # Existing
├── sprinkling.rs       # Existing
├── relations.rs        # Existing
├── causal_set.rs       # Existing + new methods
└── geometry/           # NEW - Phase 2
    ├── mod.rs          # Module exports
    ├── dimension.rs    # Myrheim-Meyer estimator
    ├── chains.rs       # Longest chain / proper time
    └── intervals.rs    # Interval counting for BD action
```

## Core Types

### DimensionEstimate

```rust
/// Result of dimension estimation from ordering fraction
#[derive(Debug, Clone, Copy)]
pub struct DimensionEstimate {
    /// Estimated spacetime dimension (NOT spatial!)
    pub dimension: f64,
    /// Whether extrapolation beyond Myrheim's Table I (d > 4) was needed
    pub extrapolated: bool,
}

impl DimensionEstimate {
    pub fn reliable(d: f64) -> Self {
        Self { dimension: d, extrapolated: false }
    }

    pub fn extrapolated(d: f64) -> Self {
        Self { dimension: d, extrapolated: true }
    }

    pub fn is_reliable(&self) -> bool {
        !self.extrapolated
    }
}

impl From<DimensionEstimate> for f64 {
    fn from(est: DimensionEstimate) -> f64 {
        est.dimension
    }
}
```

### IntervalCounts

```rust
/// Interval size distribution for Benincasa-Dowker action
#[derive(Debug, Clone)]
pub struct IntervalCounts {
    /// counts[k] = number of pairs (x,y) with exactly k elements strictly between
    pub counts: Vec<usize>,
    /// Total number of related pairs
    pub total_relations: usize,
}

impl IntervalCounts {
    pub fn n_k(&self, k: usize) -> usize {
        self.counts.get(k).copied().unwrap_or(0)
    }

    pub fn links(&self) -> usize {
        self.n_k(0)
    }

    pub fn ordering_fraction(&self, n_elements: usize) -> f64 {
        if n_elements < 2 { return 1.0; }
        2.0 * self.total_relations as f64 / (n_elements * (n_elements - 1)) as f64
    }
}
```

### DimensionStatistics

```rust
/// Statistical summary of dimension estimation
#[derive(Debug, Clone)]
pub struct DimensionStatistics {
    pub estimate: DimensionEstimate,
    /// Uncertainty δd ≈ d/√N
    pub uncertainty: f64,
    /// Ordering fraction used
    pub ordering_fraction: f64,
    /// Number of elements
    pub n_elements: usize,
}
```

## Algorithms

### Myrheim-Meyer Dimension Estimator

Uses log-linear interpolation through Myrheim 1978 Table I anchor points.

**CRITICAL NOTE:** The gamma function formula often cited in literature does NOT reproduce Myrheim's values. Use Table I as ground truth.

```rust
/// Myrheim 1978 Table I - PRIMARY SOURCE
const MYRHEIM_ANCHORS: [(f64, f64); 4] = [
    (1.0, 1.0),           // d=1 (trivial)
    (2.0, 0.5),           // d=2: f = 1/2
    (3.0, 8.0 / 35.0),    // d=3: f = 8/35 ≈ 0.2286
    (4.0, 0.1),           // d=4: f = 1/10
];

pub fn estimate_dimension(f: f64) -> DimensionEstimate {
    // Edge cases
    if f >= 1.0 { return DimensionEstimate::reliable(1.0); }
    if f <= 0.0 { return DimensionEstimate::extrapolated(f64::INFINITY); }
    if f < 0.1 { return extrapolate_high_dimension(f); }

    let ln_f = f.ln();

    // Find segment and interpolate in log-space
    for i in 0..MYRHEIM_ANCHORS.len() - 1 {
        let (d1, f1) = MYRHEIM_ANCHORS[i];
        let (d2, f2) = MYRHEIM_ANCHORS[i + 1];

        if f <= f1 && f >= f2 {
            let ln_f1 = f1.ln();
            let ln_f2 = f2.ln();
            let t = (ln_f1 - ln_f) / (ln_f1 - ln_f2);
            return DimensionEstimate::reliable(d1 + t * (d2 - d1));
        }
    }

    // Fallback (shouldn't reach here given f >= 0.1 check)
    extrapolate_high_dimension(f)
}

fn extrapolate_high_dimension(f: f64) -> DimensionEstimate {
    let ln_f = f.ln();
    let ln_f3 = (8.0_f64 / 35.0).ln();
    let ln_f4 = 0.1_f64.ln();
    let slope = ln_f4 - ln_f3;

    let d = 4.0 + (ln_f - ln_f4) / slope;
    DimensionEstimate::extrapolated(d)
}
```

### Longest Chain (Proper Time)

**Important:** Uses interval-local topological sort, not global time ordering.

```rust
pub fn longest_chain(&self, x: usize, y: usize) -> Vec<usize> {
    if !self.is_related(x, y) {
        return vec![];
    }

    let mut interval: Vec<usize> = self.interval(x, y);

    if interval.is_empty() {
        return vec![x, y];  // Direct link
    }

    // Sort interval by time (valid topo sort for sub-DAG)
    interval.sort_by(|&a, &b| {
        self.points[a].t().partial_cmp(&self.points[b].t()).unwrap()
    });

    // DP: longest[i] = (length, predecessor index)
    let mut longest: Vec<(usize, Option<usize>)> = vec![(1, None); interval.len()];

    // Base case: direct links from x
    for (i, &elem) in interval.iter().enumerate() {
        if self.is_linked(x, elem) {
            longest[i] = (2, None);
        }
    }

    // DP recurrence
    for i in 0..interval.len() {
        for j in i + 1..interval.len() {
            if self.is_related(interval[i], interval[j]) {
                let new_len = longest[i].0 + 1;
                if new_len > longest[j].0 {
                    longest[j] = (new_len, Some(i));
                }
            }
        }
    }

    // Find best endpoint with link to y
    let mut best_len = 1;
    let mut best_idx: Option<usize> = None;

    for (i, &elem) in interval.iter().enumerate() {
        if self.is_linked(elem, y) && longest[i].0 + 1 > best_len {
            best_len = longest[i].0 + 1;
            best_idx = Some(i);
        }
    }

    // Reconstruct path
    let mut chain = vec![y];
    let mut current = best_idx;
    while let Some(i) = current {
        chain.push(interval[i]);
        current = longest[i].1;
    }
    chain.push(x);
    chain.reverse();

    chain
}
```

**Complexity:** O(|interval|²)

### Interval Counting

```rust
pub fn interval_counts(&self) -> IntervalCounts {
    let n = self.len();
    let mut counts: Vec<usize> = Vec::new();
    let mut total = 0;

    for i in 0..n {
        for j in i + 1..n {
            if self.is_related(i, j) {
                let size = self.interval_size(i, j);
                if size >= counts.len() {
                    counts.resize(size + 1, 0);
                }
                counts[size] += 1;
                total += 1;
            }
        }
    }

    IntervalCounts { counts, total_relations: total }
}
```

**Complexity:** O(N³) worst case. GPU optimization deferred to Phase 3.

## New CausalSet Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `estimated_dimension()` | `DimensionEstimate` | Myrheim-Meyer from ordering fraction |
| `ordering_fraction_of_interval(x, y)` | `f64` | Ordering fraction for sub-causet |
| `longest_chain(x, y)` | `Vec<usize>` | Maximum antichain-free path |
| `longest_chain_length(x, y)` | `usize` | Length of longest chain |
| `proper_time_estimate(x, y)` | `f64` | Chain length (calibration in Phase 3) |
| `interval_counts()` | `IntervalCounts` | N_k distribution |
| `dimension_statistics()` | `DimensionStatistics` | Estimate with uncertainty |

## Validation Tests

### Anchor Point Recovery

```rust
#[test]
fn test_anchor_point_recovery() {
    assert!((estimate_dimension(0.5).dimension - 2.0).abs() < 1e-10);
    assert!((estimate_dimension(8.0/35.0).dimension - 3.0).abs() < 1e-10);
    assert!((estimate_dimension(0.1).dimension - 4.0).abs() < 1e-10);
}
```

### Monotonicity

```rust
#[test]
fn test_monotonicity() {
    let fs = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05];
    let ds: Vec<f64> = fs.iter().map(|&f| estimate_dimension(f).dimension).collect();

    for i in 1..ds.len() {
        assert!(ds[i] > ds[i-1]);
    }
}
```

### Statistical Recovery

```rust
#[test]
fn test_statistical_recovery() {
    for true_d in 2..=4 {
        let mut estimates = Vec::new();
        for _ in 0..20 {
            let causet = sprinkle_minkowski_diamond(true_d, 500);
            let f = causet.ordering_fraction();
            estimates.push(estimate_dimension(f).dimension);
        }
        let mean: f64 = estimates.iter().sum::<f64>() / estimates.len() as f64;

        assert!((mean - true_d as f64).abs() < 0.15);
    }
}
```

## Complexity Summary

| Operation | Complexity | Notes |
|-----------|------------|-------|
| `estimate_dimension()` | O(1) | 4 anchor points |
| `ordering_fraction_of_interval()` | O(interval²) | Count relations in sub-causet |
| `longest_chain()` | O(interval²) | DP on interval elements |
| `interval_counts()` | O(N³) | GPU candidate for Phase 3 |

## References

1. **Myrheim (1978)** — "Statistical Geometry" CERN-TH-2538 — Primary source for Table I
2. **Meyer (1988)** — MIT PhD thesis — Dimension estimator derivation
3. **Surya (2019)** — Living Reviews arXiv:1903.11544 — Comprehensive review
