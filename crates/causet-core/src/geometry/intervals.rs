//! Interval counting for Benincasa-Dowker action.
//!
//! The BD action uses the distribution of interval sizes:
//! N_k = number of pairs (x,y) where x ≺ y and exactly k elements lie between.
//!
//! For 4D spacetime: S^(4) ∝ N - N_1 + 9N_2 - 16N_3 + 8N_4

use crate::causal_set::CausalSet;
use serde::{Deserialize, Serialize};

/// Interval size distribution for Benincasa-Dowker action.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntervalCounts {
    /// counts[k] = number of pairs (x,y) with exactly k elements strictly between.
    /// counts[0] = N_0 (links), counts[1] = N_1, etc.
    pub counts: Vec<usize>,
    /// Total number of causally related pairs.
    pub total_relations: usize,
}

impl IntervalCounts {
    /// Create empty interval counts.
    pub fn empty() -> Self {
        Self {
            counts: Vec::new(),
            total_relations: 0,
        }
    }

    /// Get N_k (number of intervals with exactly k elements between).
    pub fn n_k(&self, k: usize) -> usize {
        self.counts.get(k).copied().unwrap_or(0)
    }

    /// Number of links (N_0) — pairs with nothing between.
    pub fn links(&self) -> usize {
        self.n_k(0)
    }

    /// Compute ordering fraction from interval counts.
    ///
    /// f = 2R / [N(N-1)] where R = total_relations
    pub fn ordering_fraction(&self, n_elements: usize) -> f64 {
        if n_elements < 2 {
            return 1.0;
        }
        2.0 * self.total_relations as f64 / (n_elements * (n_elements - 1)) as f64
    }

    /// Maximum interval size encountered.
    pub fn max_interval_size(&self) -> usize {
        if self.counts.is_empty() {
            0
        } else {
            self.counts.len() - 1
        }
    }

    /// Distribution as a vector of (size, count) pairs for non-zero counts.
    pub fn distribution(&self) -> Vec<(usize, usize)> {
        self.counts
            .iter()
            .enumerate()
            .filter(|(_, &count)| count > 0)
            .map(|(size, &count)| (size, count))
            .collect()
    }
}

impl<const D: usize> CausalSet<D> {
    /// Count intervals by size for Benincasa-Dowker action.
    ///
    /// N_k = number of pairs (x,y) where x ≺ y and exactly k elements
    /// lie strictly between them in the causal order.
    ///
    /// # Complexity
    ///
    /// O(N³) worst case — for each of O(N²) pairs, computing interval size
    /// takes O(N) time. For manifold-like causets, this is often better
    /// due to sparsity.
    ///
    /// # Example
    ///
    /// ```
    /// use causet_core::prelude::*;
    /// use causet_core::sprinkling::{CausalDiamond, Sprinkler};
    ///
    /// let diamond = CausalDiamond::<2>::symmetric(4.0);
    /// let mut sprinkler = Sprinkler::with_seed(42, 10.0);
    /// let result = sprinkler.sprinkle_diamond(&diamond);
    /// let causet = CausalSet::from_sprinkling(result);
    ///
    /// let counts = causet.interval_counts();
    /// println!("Links (N_0): {}", counts.links());
    /// println!("N_1: {}", counts.n_k(1));
    /// ```
    pub fn interval_counts(&self) -> IntervalCounts {
        let n = self.len();

        if n < 2 {
            return IntervalCounts::empty();
        }

        let mut counts: Vec<usize> = Vec::new();
        let mut total = 0usize;

        for i in 0..n {
            for j in i + 1..n {
                if self.is_related(i, j) {
                    let size = self.interval_size(i, j);

                    // Extend counts vector if needed
                    if size >= counts.len() {
                        counts.resize(size + 1, 0);
                    }

                    counts[size] += 1;
                    total += 1;
                }
            }
        }

        IntervalCounts {
            counts,
            total_relations: total,
        }
    }

    /// Compute ordering fraction for a sub-interval [x, y].
    ///
    /// This computes the ordering fraction restricted to elements
    /// in the causal interval between x and y.
    ///
    /// # Arguments
    ///
    /// * `x` - Lower bound (must be in causal past of y)
    /// * `y` - Upper bound
    ///
    /// # Returns
    ///
    /// Ordering fraction for the sub-causet, or 0.0 if x ⊀ y.
    pub fn ordering_fraction_of_interval(&self, x: usize, y: usize) -> f64 {
        if !self.is_related(x, y) {
            return 0.0;
        }

        // Get all elements in the interval (including x and y)
        let mut elements: Vec<usize> = self.interval(x, y);
        elements.push(x);
        elements.push(y);

        let n = elements.len();
        if n < 2 {
            return 1.0;
        }

        // Count relations within this subset
        let mut relations = 0usize;
        for i in 0..elements.len() {
            for j in i + 1..elements.len() {
                if self.is_related(elements[i], elements[j])
                    || self.is_related(elements[j], elements[i])
                {
                    relations += 1;
                }
            }
        }

        2.0 * relations as f64 / (n * (n - 1)) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point2D;

    #[test]
    fn test_interval_counts_chain() {
        // Total order: 0 ≺ 1 ≺ 2 ≺ 3 ≺ 4
        let points: Vec<Point2D> = (0..5)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);
        let counts = causet.interval_counts();

        // In a chain of 5:
        // Links (N_0): (0,1), (1,2), (2,3), (3,4) = 4
        // N_1: (0,2), (1,3), (2,4) = 3  (one element between)
        // N_2: (0,3), (1,4) = 2
        // N_3: (0,4) = 1

        assert_eq!(counts.links(), 4, "Expected 4 links");
        assert_eq!(counts.n_k(1), 3, "Expected 3 pairs with 1 between");
        assert_eq!(counts.n_k(2), 2, "Expected 2 pairs with 2 between");
        assert_eq!(counts.n_k(3), 1, "Expected 1 pair with 3 between");

        // Total relations: 4 + 3 + 2 + 1 = 10
        assert_eq!(counts.total_relations, 10);
    }

    #[test]
    fn test_interval_counts_antichain() {
        // All spacelike separated (antichain)
        let points = vec![
            Point2D::new_2d(0.0, -2.0, 0),
            Point2D::new_2d(0.0, 0.0, 1),
            Point2D::new_2d(0.0, 2.0, 2),
        ];

        let causet = CausalSet::from_points(points);
        let counts = causet.interval_counts();

        assert_eq!(counts.total_relations, 0);
        assert!(counts.counts.is_empty());
    }

    #[test]
    fn test_interval_counts_diamond() {
        // Diamond pattern: 0 → {1,2} → 3
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0),
            Point2D::new_2d(1.0, 0.3, 1),
            Point2D::new_2d(1.0, -0.3, 2),
            Point2D::new_2d(2.0, 0.0, 3),
        ];

        let causet = CausalSet::from_points(points);
        let counts = causet.interval_counts();

        // Relations:
        // 0≺1, 0≺2, 0≺3 (0 is in past of all)
        // 1≺3, 2≺3 (1,2 are in past of 3)
        // 1 and 2 are spacelike
        // Links: 0→1, 0→2, 1→3, 2→3 = 4
        // N_1: 0→3 (has 1 or 2 between, but which counts?)
        // Actually interval(0,3) = {1, 2}, so N_2: 1 pair with 2 between

        // Total: 5 relations
        assert_eq!(counts.total_relations, 5);
        assert_eq!(counts.links(), 4); // Direct links
    }

    #[test]
    fn test_ordering_fraction_of_interval() {
        // Chain of 5
        let points: Vec<Point2D> = (0..5)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);

        // Interval [0, 4] includes all elements, should be close to global f
        let f_interval = causet.ordering_fraction_of_interval(0, 4);
        let _f_global = causet.ordering_fraction();

        // In a total order, f should be 1.0
        assert!(
            (f_interval - 1.0).abs() < 0.01,
            "Expected ~1.0, got {}",
            f_interval
        );
    }

    #[test]
    fn test_ordering_fraction_unrelated() {
        let points = vec![
            Point2D::new_2d(0.0, -2.0, 0),
            Point2D::new_2d(0.0, 2.0, 1),
        ];

        let causet = CausalSet::from_points(points);

        // Spacelike separated points
        assert_eq!(causet.ordering_fraction_of_interval(0, 1), 0.0);
    }

    #[test]
    fn test_interval_counts_consistency() {
        // Verify that total_relations matches sum of counts
        let points: Vec<Point2D> = (0..10)
            .map(|i| Point2D::new_2d(i as f64, (i as f64 * 0.1).sin(), i))
            .collect();

        let causet = CausalSet::from_points(points);
        let counts = causet.interval_counts();

        let sum: usize = counts.counts.iter().sum();
        assert_eq!(
            sum, counts.total_relations,
            "Sum of counts should equal total_relations"
        );
    }

    #[test]
    fn test_distribution() {
        let points: Vec<Point2D> = (0..5)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);
        let counts = causet.interval_counts();

        let dist = counts.distribution();
        assert!(!dist.is_empty());

        // All entries should have non-zero counts
        for (_, count) in &dist {
            assert!(*count > 0);
        }
    }
}
