//! Myrheim-Meyer dimension estimator.
//!
//! Extracts spacetime dimension from causal structure using the ordering fraction.
//! Uses Myrheim 1978 Table I as ground truth with log-linear interpolation.
//!
//! # Important Note
//!
//! The gamma function formula often cited in literature:
//! `f(d) = Γ(d+1) · Γ(d/2) / [4 · Γ(3d/2)]`
//! does NOT reproduce Myrheim's original values. This implementation uses
//! the exact values from Table I of CERN-TH-2538 (1978).

use serde::{Deserialize, Serialize};

/// Myrheim 1978 Table I - PRIMARY SOURCE
/// (spacetime_dimension, ordering_fraction)
///
/// These are the AUTHORITATIVE values. Any other values in literature
/// that differ from these are propagated errors.
const MYRHEIM_ANCHORS: [(f64, f64); 4] = [
    (1.0, 1.0),        // d=1 spacetime (trivial, totally ordered)
    (2.0, 0.5),        // d=2 spacetime (1+1D): f = 1/2
    (3.0, 8.0 / 35.0), // d=3 spacetime (2+1D): f = 8/35 ≈ 0.2286
    (4.0, 0.1),        // d=4 spacetime (3+1D): f = 1/10
];

/// Result of dimension estimation from ordering fraction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DimensionEstimate {
    /// Estimated spacetime dimension (NOT spatial dimensions!)
    pub dimension: f64,
    /// Whether extrapolation beyond Myrheim's Table I (d > 4) was needed
    pub extrapolated: bool,
}

impl DimensionEstimate {
    /// Create a reliable estimate (within validated range d ∈ [1, 4]).
    pub fn reliable(d: f64) -> Self {
        Self {
            dimension: d,
            extrapolated: false,
        }
    }

    /// Create an extrapolated estimate (d > 4, less reliable).
    pub fn extrapolated(d: f64) -> Self {
        Self {
            dimension: d,
            extrapolated: true,
        }
    }

    /// Is this estimate within the validated range?
    pub fn is_reliable(&self) -> bool {
        !self.extrapolated
    }
}

impl From<DimensionEstimate> for f64 {
    fn from(est: DimensionEstimate) -> f64 {
        est.dimension
    }
}

/// Statistical summary of dimension estimation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DimensionStatistics {
    /// The dimension estimate
    pub estimate: DimensionEstimate,
    /// Uncertainty δd ≈ d/√N (statistical sampling error)
    pub uncertainty: f64,
    /// Ordering fraction used for estimation
    pub ordering_fraction: f64,
    /// Number of elements in the causal set
    pub n_elements: usize,
}

impl DimensionStatistics {
    /// Compute dimension statistics from ordering fraction and element count.
    pub fn compute(ordering_fraction: f64, n_elements: usize) -> Self {
        let estimate = estimate_dimension(ordering_fraction);

        // Statistical uncertainty: δd/d ≈ 1/√N
        // This comes from the variance in the ordering fraction estimate
        let uncertainty = if n_elements > 1 {
            estimate.dimension / (n_elements as f64).sqrt()
        } else {
            f64::INFINITY
        };

        Self {
            estimate,
            uncertainty,
            ordering_fraction,
            n_elements,
        }
    }

    /// 95% confidence interval for dimension estimate.
    pub fn confidence_interval_95(&self) -> (f64, f64) {
        let delta = 1.96 * self.uncertainty;
        (self.estimate.dimension - delta, self.estimate.dimension + delta)
    }
}

/// Estimate spacetime dimension from observed ordering fraction.
///
/// Uses log-linear interpolation through Myrheim 1978 Table I anchor points.
/// The ordering fraction f = 2R / [N(N-1)] where R is the number of
/// causally related pairs.
///
/// # Arguments
///
/// * `f` - Observed ordering fraction (should be in range (0, 1])
///
/// # Returns
///
/// `DimensionEstimate` with the estimated dimension and reliability flag.
///
/// # Example
///
/// ```
/// use causet_core::geometry::estimate_dimension;
///
/// // 2D spacetime has ordering fraction 0.5
/// let est = estimate_dimension(0.5);
/// assert!((est.dimension - 2.0).abs() < 1e-10);
/// assert!(est.is_reliable());
///
/// // Ordering fraction below 0.1 requires extrapolation
/// let est = estimate_dimension(0.05);
/// assert!(!est.is_reliable());
/// ```
pub fn estimate_dimension(f: f64) -> DimensionEstimate {
    // Edge cases
    if f >= 1.0 {
        return DimensionEstimate::reliable(1.0);
    }
    if f <= 0.0 {
        return DimensionEstimate::extrapolated(f64::INFINITY);
    }
    if f < 0.1 {
        // Below d=4 anchor point, need extrapolation
        return extrapolate_high_dimension(f);
    }

    let ln_f = f.ln();

    // Find the segment containing f and interpolate in log-space
    for i in 0..MYRHEIM_ANCHORS.len() - 1 {
        let (d1, f1) = MYRHEIM_ANCHORS[i];
        let (d2, f2) = MYRHEIM_ANCHORS[i + 1];

        if f <= f1 && f >= f2 {
            let ln_f1 = f1.ln();
            let ln_f2 = f2.ln();
            // Linear interpolation in log(f) space
            let t = (ln_f1 - ln_f) / (ln_f1 - ln_f2);
            return DimensionEstimate::reliable(d1 + t * (d2 - d1));
        }
    }

    // Fallback (shouldn't reach here given f >= 0.1 check above)
    extrapolate_high_dimension(f)
}

/// Extrapolate dimension for f < 0.1 (d > 4).
///
/// Continues the log-linear trend from the d=3→4 segment.
fn extrapolate_high_dimension(f: f64) -> DimensionEstimate {
    let ln_f = f.ln();
    let ln_f3 = (8.0_f64 / 35.0).ln(); // ≈ -1.475
    let ln_f4 = 0.1_f64.ln(); // ≈ -2.303
    let slope = ln_f4 - ln_f3; // change in ln(f) per unit dimension

    let d = 4.0 + (ln_f - ln_f4) / slope;
    DimensionEstimate::extrapolated(d)
}

/// Get the theoretical ordering fraction for a given spacetime dimension.
///
/// Uses log-linear interpolation between Myrheim's anchor points.
/// This is the inverse of `estimate_dimension`.
pub fn theoretical_ordering_fraction(d: f64) -> f64 {
    if d <= 1.0 {
        return 1.0;
    }
    if d >= 4.0 {
        // Extrapolate using log-linear trend
        let ln_f3 = (8.0_f64 / 35.0).ln();
        let ln_f4 = 0.1_f64.ln();
        let slope = ln_f4 - ln_f3;
        return (ln_f4 + slope * (d - 4.0)).exp();
    }

    // Find segment and interpolate
    for i in 0..MYRHEIM_ANCHORS.len() - 1 {
        let (d1, f1) = MYRHEIM_ANCHORS[i];
        let (d2, f2) = MYRHEIM_ANCHORS[i + 1];

        if d >= d1 && d <= d2 {
            let t = (d - d1) / (d2 - d1);
            let ln_f1 = f1.ln();
            let ln_f2 = f2.ln();
            return (ln_f1 + t * (ln_f2 - ln_f1)).exp();
        }
    }

    // Shouldn't reach here
    0.1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_anchor_point_recovery() {
        // Should recover exact dimensions at anchor points
        assert!((estimate_dimension(1.0).dimension - 1.0).abs() < 1e-10);
        assert!((estimate_dimension(0.5).dimension - 2.0).abs() < 1e-10);
        assert!((estimate_dimension(8.0 / 35.0).dimension - 3.0).abs() < 1e-10);
        assert!((estimate_dimension(0.1).dimension - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_anchor_points_reliable() {
        assert!(estimate_dimension(1.0).is_reliable());
        assert!(estimate_dimension(0.5).is_reliable());
        assert!(estimate_dimension(8.0 / 35.0).is_reliable());
        assert!(estimate_dimension(0.1).is_reliable());
    }

    #[test]
    fn test_extrapolation_flagged() {
        // Below f=0.1 should be flagged as extrapolated
        assert!(!estimate_dimension(0.09).is_reliable());
        assert!(!estimate_dimension(0.05).is_reliable());
        assert!(!estimate_dimension(0.01).is_reliable());
    }

    #[test]
    fn test_monotonicity() {
        // Higher ordering fraction → lower dimension
        let fs = [0.9, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01];
        let ds: Vec<f64> = fs.iter().map(|&f| estimate_dimension(f).dimension).collect();

        for i in 1..ds.len() {
            assert!(
                ds[i] > ds[i - 1],
                "Monotonicity violated: d({})={}, d({})={}",
                fs[i - 1],
                ds[i - 1],
                fs[i],
                ds[i]
            );
        }
    }

    #[test]
    fn test_interpolation_midpoints() {
        // Test that midpoint ordering fractions give midpoint dimensions (roughly)
        // Between d=2 and d=3
        let f_mid = ((0.5_f64).ln() + (8.0 / 35.0_f64).ln()) / 2.0;
        let f_mid = f_mid.exp();
        let d_mid = estimate_dimension(f_mid).dimension;
        assert!(
            (d_mid - 2.5).abs() < 0.1,
            "Expected ~2.5, got {}",
            d_mid
        );
    }

    #[test]
    fn test_theoretical_ordering_fraction_inverse() {
        // theoretical_ordering_fraction should be inverse of estimate_dimension
        for d in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0] {
            let f = theoretical_ordering_fraction(d);
            let d_recovered = estimate_dimension(f).dimension;
            assert!(
                (d_recovered - d).abs() < 1e-6,
                "d={}: f={}, d_recovered={}",
                d,
                f,
                d_recovered
            );
        }
    }

    #[test]
    fn test_edge_cases() {
        assert_eq!(estimate_dimension(1.0).dimension, 1.0);
        assert_eq!(estimate_dimension(1.5).dimension, 1.0); // Clamped
        assert!(estimate_dimension(0.0).dimension.is_infinite());
        assert!(estimate_dimension(-0.1).dimension.is_infinite());
    }

    #[test]
    fn test_dimension_statistics() {
        let stats = DimensionStatistics::compute(0.5, 1000);
        assert!((stats.estimate.dimension - 2.0).abs() < 1e-10);
        assert!(stats.uncertainty < 0.1); // ~2/√1000 ≈ 0.063
        assert!(stats.estimate.is_reliable());
    }

    #[test]
    fn test_confidence_interval() {
        let stats = DimensionStatistics::compute(0.5, 1000);
        let (lo, hi) = stats.confidence_interval_95();
        assert!(lo < 2.0);
        assert!(hi > 2.0);
        assert!(hi - lo < 0.3); // Reasonable interval width
    }
}
