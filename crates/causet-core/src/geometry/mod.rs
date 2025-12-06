//! Geometric observables for causal sets.
//!
//! This module implements Phase 2 of the CAUSET project: extracting spacetime
//! geometry from pure causal structure.
//!
//! # Key Components
//!
//! - **Dimension estimation**: Myrheim-Meyer estimator recovers spacetime
//!   dimension from the ordering fraction.
//! - **Longest chains**: Proper time estimation via geodesic-like paths.
//! - **Interval counting**: N_k distribution for Benincasa-Dowker action.
//!
//! # Example
//!
//! ```
//! use causet_core::prelude::*;
//! use causet_core::geometry::{estimate_dimension, DimensionStatistics};
//! use causet_core::sprinkling::{CausalDiamond, Sprinkler};
//!
//! // Sprinkle into 2D Minkowski spacetime
//! let diamond = CausalDiamond::<2>::symmetric(4.0);
//! let mut sprinkler = Sprinkler::with_seed(42, 50.0);
//! let result = sprinkler.sprinkle_diamond(&diamond);
//! let causet = CausalSet::from_sprinkling(result);
//!
//! // Estimate dimension from ordering fraction
//! let f = causet.ordering_fraction();
//! let dim = estimate_dimension(f);
//! println!("Ordering fraction: {:.4}", f);
//! println!("Estimated dimension: {:.2}", dim.dimension);
//!
//! // Get full statistics
//! let stats = causet.dimension_statistics();
//! println!("Dimension: {:.2} ± {:.2}", stats.estimate.dimension, stats.uncertainty);
//! ```

mod chains;
mod dimension;
mod intervals;

// Re-export main types and functions
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
}
