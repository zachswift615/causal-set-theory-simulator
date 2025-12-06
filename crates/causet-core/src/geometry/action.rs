//! Benincasa-Dowker action for causal sets.
//!
//! The BD action is a discrete analog of the Einstein-Hilbert action that
//! recovers general relativity in the continuum limit. It encodes scalar
//! curvature through interval counting.
//!
//! # Formulas (verified from Benincasa-Dowker 2010, arXiv:1001.2725v4)
//!
//! **4D action:**
//! ```text
//! S^(4) = N - N₁ + 9N₂ - 16N₃ + 8N₄
//! ```
//!
//! **2D action:**
//! ```text
//! S^(2) = N - 2N₁ + 4N₂ - 2N₃
//! ```
//!
//! Where N_k counts (k+1)-element **inclusive** order intervals:
//! - N₁ = links (pairs with 0 elements strictly between)
//! - N₂ = pairs with 1 element strictly between
//! - N₃ = pairs with 2 elements strictly between
//! - N₄ = pairs with 3 elements strictly between
//!
//! # Physical Interpretation
//!
//! For flat Minkowski spacetime (zero curvature), the action should
//! fluctuate around zero with magnitude O(√N).
//!
//! For curved spacetimes, the action encodes the integrated scalar curvature.

use super::IntervalCounts;
use serde::{Deserialize, Serialize};

/// Benincasa-Dowker action coefficients for different dimensions.
///
/// The action formula is: S = N + Σ c_k · N_k
/// where c_k are the coefficients and N_k count interval sizes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BDCoefficients {
    /// Spacetime dimension
    pub dimension: usize,
    /// Coefficients for N_1, N_2, N_3, ... (N_k = pairs with k-1 elements between)
    /// Index 0 = coefficient for links (N_1)
    pub coefficients: &'static [f64],
}

/// 2D (1+1) Benincasa-Dowker coefficients
/// S^(2) = N - 2N₁ + 4N₂ - 2N₃
pub const BD_COEFFICIENTS_2D: BDCoefficients = BDCoefficients {
    dimension: 2,
    coefficients: &[-2.0, 4.0, -2.0],
};

/// 4D (3+1) Benincasa-Dowker coefficients
/// S^(4) = N - N₁ + 9N₂ - 16N₃ + 8N₄
pub const BD_COEFFICIENTS_4D: BDCoefficients = BDCoefficients {
    dimension: 4,
    coefficients: &[-1.0, 9.0, -16.0, 8.0],
};

/// Result of BD action computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BDActionResult {
    /// The computed action value
    pub action: f64,
    /// Number of elements in the causal set
    pub n_elements: usize,
    /// Spacetime dimension used
    pub dimension: usize,
    /// Action per element (useful for comparing different sizes)
    pub action_per_element: f64,
    /// Expected fluctuation magnitude for flat spacetime: O(√N)
    pub expected_fluctuation: f64,
    /// Interval counts used in computation
    pub interval_counts: IntervalCounts,
}

impl BDActionResult {
    /// Is the action consistent with flat spacetime (zero curvature)?
    ///
    /// Returns true if |S| < k * √N where k is typically 2-3 sigma.
    pub fn is_consistent_with_flat(&self, sigma: f64) -> bool {
        self.action.abs() < sigma * self.expected_fluctuation
    }

    /// Normalized action: S / √N
    ///
    /// For flat spacetime, this should be O(1) regardless of N.
    pub fn normalized_action(&self) -> f64 {
        if self.n_elements > 0 {
            self.action / (self.n_elements as f64).sqrt()
        } else {
            0.0
        }
    }
}

/// Compute the Benincasa-Dowker action for 2D spacetime.
///
/// Formula: S^(2) = N - 2N₁ + 4N₂ - 2N₃
///
/// Verified from Benincasa-Dowker 2010 (arXiv:1001.2725v4), line 599.
pub fn bd_action_2d(counts: &IntervalCounts, n_elements: usize) -> BDActionResult {
    // N_k in BD notation = our counts[k-1]
    // N₁ = links = counts[0]
    // N₂ = 1 between = counts[1]
    // N₃ = 2 between = counts[2]
    let n1 = counts.n_k(0) as f64;
    let n2 = counts.n_k(1) as f64;
    let n3 = counts.n_k(2) as f64;

    let n = n_elements as f64;
    let action = n - 2.0 * n1 + 4.0 * n2 - 2.0 * n3;

    BDActionResult {
        action,
        n_elements,
        dimension: 2,
        action_per_element: action / n,
        expected_fluctuation: n.sqrt(),
        interval_counts: counts.clone(),
    }
}

/// Compute the Benincasa-Dowker action for 4D spacetime.
///
/// Formula: S^(4) = N - N₁ + 9N₂ - 16N₃ + 8N₄
///
/// Verified from Benincasa-Dowker 2010 (arXiv:1001.2725v4), lines 602-604.
pub fn bd_action_4d(counts: &IntervalCounts, n_elements: usize) -> BDActionResult {
    // N_k in BD notation = our counts[k-1]
    let n1 = counts.n_k(0) as f64; // links
    let n2 = counts.n_k(1) as f64; // 1 between
    let n3 = counts.n_k(2) as f64; // 2 between
    let n4 = counts.n_k(3) as f64; // 3 between

    let n = n_elements as f64;
    let action = n - n1 + 9.0 * n2 - 16.0 * n3 + 8.0 * n4;

    BDActionResult {
        action,
        n_elements,
        dimension: 4,
        action_per_element: action / n,
        expected_fluctuation: n.sqrt(),
        interval_counts: counts.clone(),
    }
}

/// Compute the BD action using generic coefficients.
///
/// For dimensions other than 2D and 4D, you can provide custom coefficients.
pub fn bd_action_generic(
    counts: &IntervalCounts,
    n_elements: usize,
    coefficients: &BDCoefficients,
) -> BDActionResult {
    let n = n_elements as f64;

    // Start with N
    let mut action = n;

    // Add coefficient * N_k for each k
    for (k, &coeff) in coefficients.coefficients.iter().enumerate() {
        action += coeff * counts.n_k(k) as f64;
    }

    BDActionResult {
        action,
        n_elements,
        dimension: coefficients.dimension,
        action_per_element: action / n,
        expected_fluctuation: n.sqrt(),
        interval_counts: counts.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bd_action_2d_formula() {
        // Test with known counts
        let counts = IntervalCounts {
            counts: vec![100, 50, 20], // N₁=100, N₂=50, N₃=20
            total_relations: 170,
        };

        let result = bd_action_2d(&counts, 500);

        // S = 500 - 2*100 + 4*50 - 2*20 = 500 - 200 + 200 - 40 = 460
        assert_eq!(result.action, 460.0);
        assert_eq!(result.dimension, 2);
    }

    #[test]
    fn test_bd_action_4d_formula() {
        // Test with known counts
        let counts = IntervalCounts {
            counts: vec![100, 50, 20, 10], // N₁=100, N₂=50, N₃=20, N₄=10
            total_relations: 180,
        };

        let result = bd_action_4d(&counts, 500);

        // S = 500 - 100 + 9*50 - 16*20 + 8*10
        //   = 500 - 100 + 450 - 320 + 80 = 610
        assert_eq!(result.action, 610.0);
        assert_eq!(result.dimension, 4);
    }

    #[test]
    fn test_bd_action_generic_matches_specific() {
        let counts = IntervalCounts {
            counts: vec![100, 50, 20, 10],
            total_relations: 180,
        };

        let specific = bd_action_4d(&counts, 500);
        let generic = bd_action_generic(&counts, 500, &BD_COEFFICIENTS_4D);

        assert!((specific.action - generic.action).abs() < 1e-10);
    }

    #[test]
    fn test_normalized_action() {
        let counts = IntervalCounts {
            counts: vec![100, 50, 20],
            total_relations: 170,
        };

        let result = bd_action_2d(&counts, 500);

        // normalized = S / √N = 460 / √500 ≈ 20.6
        let expected = 460.0 / 500.0_f64.sqrt();
        assert!((result.normalized_action() - expected).abs() < 0.01);
    }

    #[test]
    fn test_flat_spacetime_consistency() {
        // For flat spacetime, action should be small relative to √N
        let counts = IntervalCounts {
            counts: vec![100, 50, 20],
            total_relations: 170,
        };

        let result = bd_action_2d(&counts, 500);

        // This particular example won't be "flat" but test the method
        assert!(result.expected_fluctuation > 0.0);
        // is_consistent_with_flat checks |S| < sigma * √N
    }

    #[test]
    fn test_empty_causet() {
        let counts = IntervalCounts::empty();
        let result = bd_action_4d(&counts, 0);

        assert!(result.action.is_nan() || result.action == 0.0);
    }
}
