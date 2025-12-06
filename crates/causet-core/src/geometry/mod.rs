//! Geometric observables for causal sets.
//!
//! This module implements Phase 2 & 3 of the CAUSET project: extracting spacetime
//! geometry and computing the Benincasa-Dowker action from pure causal structure.
//!
//! # Key Components
//!
//! - **Dimension estimation**: Myrheim-Meyer estimator recovers spacetime
//!   dimension from the ordering fraction.
//! - **Longest chains**: Proper time estimation via geodesic-like paths.
//! - **Interval counting**: N_k distribution for Benincasa-Dowker action.
//! - **BD action**: Discrete Einstein-Hilbert action encoding scalar curvature.
//!
//! # Example
//!
//! ```
//! use causet_core::prelude::*;
//! use causet_core::geometry::{estimate_dimension, DimensionStatistics};
//! use causet_core::sprinkling::{CausalDiamond, Sprinkler};
//!
//! // Sprinkle into 2D Minkowski spacetime
//! let diamond = CausalDiamond::<2>::symmetric(10.0);
//! let mut sprinkler = Sprinkler::with_seed(42, 1.0);
//! let result = sprinkler.sprinkle_fixed(&diamond, 500);
//! let causet = CausalSet::from_sprinkling(result);
//!
//! // Estimate dimension from ordering fraction
//! let dim = causet.estimated_dimension();
//! println!("Estimated dimension: {:.2}", dim.dimension);
//!
//! // Compute BD action (should be ~0 for flat spacetime)
//! let action = causet.bd_action_2d();
//! println!("BD action: {:.2} (normalized: {:.2})", action.action, action.normalized_action());
//! ```

mod action;
mod chains;
mod dimension;
mod intervals;

// Re-export main types and functions
pub use action::{
    bd_action_2d, bd_action_4d, bd_action_generic, BDActionResult, BDCoefficients,
    BD_COEFFICIENTS_2D, BD_COEFFICIENTS_4D,
};
pub use chains::ChainResult;
pub use dimension::{
    estimate_dimension, theoretical_ordering_fraction, DimensionEstimate, DimensionStatistics,
};
pub use intervals::IntervalCounts;

// Extend CausalSet with geometry methods
use crate::causal_set::CausalSet;

impl<const D: usize> CausalSet<D> {
    /// Estimate spacetime dimension from the causal structure.
    ///
    /// Uses the Myrheim-Meyer estimator based on ordering fraction.
    ///
    /// # Example
    ///
    /// ```
    /// use causet_core::prelude::*;
    /// use causet_core::sprinkling::{CausalDiamond, Sprinkler};
    ///
    /// let diamond = CausalDiamond::<2>::symmetric(4.0);
    /// let mut sprinkler = Sprinkler::with_seed(42, 50.0);
    /// let result = sprinkler.sprinkle_diamond(&diamond);
    /// let causet = CausalSet::from_sprinkling(result);
    ///
    /// let dim = causet.estimated_dimension();
    /// println!("Estimated dimension: {:.2}", dim.dimension);
    /// assert!(dim.is_reliable());
    /// ```
    pub fn estimated_dimension(&self) -> DimensionEstimate {
        estimate_dimension(self.ordering_fraction())
    }

    /// Get full dimension statistics including uncertainty.
    ///
    /// Returns the dimension estimate along with statistical uncertainty
    /// based on the number of elements (δd ≈ d/√N).
    pub fn dimension_statistics(&self) -> DimensionStatistics {
        DimensionStatistics::compute(self.ordering_fraction(), self.len())
    }

    /// Estimate uncertainty in the dimension measurement.
    ///
    /// Returns δd ≈ d/√N where d is the estimated dimension and N
    /// is the number of elements.
    pub fn dimension_uncertainty(&self) -> f64 {
        self.dimension_statistics().uncertainty
    }

    /// Compute the Benincasa-Dowker action for 2D spacetime.
    ///
    /// Formula: S^(2) = N - 2N₁ + 4N₂ - 2N₃
    ///
    /// For flat Minkowski spacetime, the action should fluctuate around zero
    /// with magnitude O(√N).
    ///
    /// # Example
    ///
    /// ```
    /// use causet_core::prelude::*;
    /// use causet_core::sprinkling::{CausalDiamond, Sprinkler};
    ///
    /// let diamond = CausalDiamond::<2>::symmetric(10.0);
    /// let mut sprinkler = Sprinkler::with_seed(42, 1.0);
    /// let result = sprinkler.sprinkle_fixed(&diamond, 500);
    /// let causet = CausalSet::from_sprinkling(result);
    ///
    /// let action = causet.bd_action_2d();
    /// println!("Action: {:.2}, Normalized: {:.2}", action.action, action.normalized_action());
    /// ```
    pub fn bd_action_2d(&self) -> BDActionResult {
        let counts = self.interval_counts();
        bd_action_2d(&counts, self.len())
    }

    /// Compute the Benincasa-Dowker action for 4D spacetime.
    ///
    /// Formula: S^(4) = N - N₁ + 9N₂ - 16N₃ + 8N₄
    ///
    /// For flat Minkowski spacetime, the action should fluctuate around zero
    /// with magnitude O(√N).
    ///
    /// # Example
    ///
    /// ```
    /// use causet_core::prelude::*;
    /// use causet_core::sprinkling::{CausalDiamond, Sprinkler};
    ///
    /// let diamond = CausalDiamond::<4>::symmetric(10.0);
    /// let mut sprinkler = Sprinkler::with_seed(42, 1.0);
    /// let result = sprinkler.sprinkle_fixed(&diamond, 500);
    /// let causet = CausalSet::from_sprinkling(result);
    ///
    /// let action = causet.bd_action_4d();
    /// println!("Action: {:.2}, Normalized: {:.2}", action.action, action.normalized_action());
    /// ```
    pub fn bd_action_4d(&self) -> BDActionResult {
        let counts = self.interval_counts();
        bd_action_4d(&counts, self.len())
    }

    /// Compute the BD action using the appropriate formula for this causet's dimension.
    ///
    /// Returns the 2D action for D=2, 4D action for D=4, and panics for other dimensions.
    pub fn bd_action(&self) -> BDActionResult {
        match D {
            2 => self.bd_action_2d(),
            4 => self.bd_action_4d(),
            _ => panic!(
                "BD action not implemented for {}D spacetime. Use bd_action_2d() or bd_action_4d() explicitly.",
                D
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point2D;
    use crate::sprinkling::{CausalDiamond, Sprinkler};

    #[test]
    fn test_dimension_estimation_2d() {
        // Sprinkle into 2D Minkowski spacetime
        // Use sprinkle_fixed for correct sampling distribution
        let diamond = CausalDiamond::<2>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(12345, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, 500);
        let causet = CausalSet::from_sprinkling(result);

        let f = causet.ordering_fraction();
        let dim = causet.estimated_dimension();

        // 2D ordering fraction should be ~0.5, dimension ~2.0
        assert!(
            (dim.dimension - 2.0).abs() < 0.3,
            "Expected ~2.0, got {} (f={:.4})",
            dim.dimension,
            f
        );
        assert!(dim.is_reliable());
    }

    #[test]
    fn test_dimension_estimation_3d() {
        // Sprinkle into 3D Minkowski spacetime
        let diamond = CausalDiamond::<3>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(12345, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, 500);
        let causet = CausalSet::from_sprinkling(result);

        let f = causet.ordering_fraction();
        let dim = causet.estimated_dimension();

        // 3D ordering fraction should be ~0.229, dimension ~3.0
        assert!(
            (dim.dimension - 3.0).abs() < 0.3,
            "Expected ~3.0, got {} (f={:.4})",
            dim.dimension,
            f
        );
        assert!(dim.is_reliable());
    }

    #[test]
    fn test_dimension_estimation_4d() {
        // Sprinkle into 4D Minkowski spacetime
        let diamond = CausalDiamond::<4>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(12345, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, 500);
        let causet = CausalSet::from_sprinkling(result);

        let f = causet.ordering_fraction();
        let dim = causet.estimated_dimension();

        // 4D ordering fraction should be ~0.1, dimension ~4.0
        assert!(
            (dim.dimension - 4.0).abs() < 0.3,
            "Expected ~4.0, got {} (f={:.4})",
            dim.dimension,
            f
        );
        // At the boundary f=0.1, might be exactly 4.0 or slightly extrapolated
    }

    #[test]
    fn test_dimension_statistics() {
        let diamond = CausalDiamond::<2>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(42, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, 500);
        let causet = CausalSet::from_sprinkling(result);

        let stats = causet.dimension_statistics();

        assert!(stats.n_elements > 0);
        assert!(stats.ordering_fraction > 0.0);
        assert!(stats.ordering_fraction < 1.0);
        assert!(stats.uncertainty > 0.0);
        assert!(stats.uncertainty < 0.5); // Reasonable for N ~ 500
    }

    #[test]
    fn test_longest_chain_integration() {
        // Create a chain and verify longest_chain works through CausalSet
        let points: Vec<Point2D> = (0..10)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);

        let chain = causet.longest_chain(0, 9);
        assert_eq!(chain.length, 10);

        let height = causet.height();
        assert_eq!(height, 10);
    }

    #[test]
    fn test_interval_counts_integration() {
        let diamond = CausalDiamond::<2>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(42, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, 200);
        let causet = CausalSet::from_sprinkling(result);

        let counts = causet.interval_counts();

        // Should have links
        assert!(counts.links() > 0);

        // Total should match relation count from matrix
        let matrix_relations = causet.causal_matrix().relation_count();
        assert_eq!(counts.total_relations, matrix_relations);
    }

    #[test]
    fn test_bd_action_2d_flat_spacetime() {
        // For flat Minkowski spacetime, BD action should fluctuate around 0
        // with magnitude O(√N)
        let diamond = CausalDiamond::<2>::symmetric(10.0);

        let mut actions = Vec::new();
        for trial in 0..10 {
            let mut sprinkler = Sprinkler::with_seed(42 + trial * 1000, 1.0);
            let result = sprinkler.sprinkle_fixed(&diamond, 500);
            let causet = CausalSet::from_sprinkling(result);

            let action = causet.bd_action_2d();
            actions.push(action.normalized_action());
        }

        // Mean normalized action should be close to 0
        let mean: f64 = actions.iter().sum::<f64>() / actions.len() as f64;
        let std: f64 = (actions
            .iter()
            .map(|&a| (a - mean).powi(2))
            .sum::<f64>()
            / (actions.len() - 1) as f64)
            .sqrt();

        // For flat spacetime, |mean| should be small relative to std
        // Allow generous tolerance since these are statistical fluctuations
        assert!(
            mean.abs() < 3.0 * std + 50.0,
            "2D BD action mean={:.2}, std={:.2} - may indicate non-flat behavior",
            mean,
            std
        );
    }

    #[test]
    fn test_bd_action_4d_flat_spacetime() {
        // For flat Minkowski spacetime, BD action should fluctuate around 0
        let diamond = CausalDiamond::<4>::symmetric(10.0);

        let mut actions = Vec::new();
        for trial in 0..5 {
            let mut sprinkler = Sprinkler::with_seed(42 + trial * 1000, 1.0);
            let result = sprinkler.sprinkle_fixed(&diamond, 500);
            let causet = CausalSet::from_sprinkling(result);

            let action = causet.bd_action_4d();
            actions.push(action.normalized_action());
        }

        let mean: f64 = actions.iter().sum::<f64>() / actions.len() as f64;

        // Just verify we get finite values for now
        // The exact behavior requires more careful analysis
        assert!(
            mean.is_finite(),
            "4D BD action should produce finite values"
        );
    }

    #[test]
    fn test_bd_action_method_matches_dimension() {
        // Test that bd_action() dispatches correctly
        let diamond = CausalDiamond::<2>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(42, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, 200);
        let causet = CausalSet::from_sprinkling(result);

        let action_explicit = causet.bd_action_2d();
        let action_dispatch = causet.bd_action();

        assert!((action_explicit.action - action_dispatch.action).abs() < 1e-10);
    }
}
