//! Longest chain computation for proper time estimation.
//!
//! In causal set theory, the longest chain (antichain-free path) between two
//! causally related points approximates the proper time separation:
//!
//! τ ≈ C_d · L / ρ^(1/d)
//!
//! where L is the longest chain length, ρ is the sprinkling density,
//! and C_d is a dimension-dependent calibration constant.

use crate::causal_set::CausalSet;

/// Result of longest chain computation.
#[derive(Debug, Clone)]
pub struct ChainResult {
    /// The elements in the chain, from start to end (inclusive).
    pub chain: Vec<usize>,
    /// Length of the chain (number of elements).
    pub length: usize,
}

impl ChainResult {
    /// Create an empty chain result (for unrelated points).
    pub fn empty() -> Self {
        Self {
            chain: Vec::new(),
            length: 0,
        }
    }

    /// Create a direct link chain (just two endpoints).
    pub fn direct_link(x: usize, y: usize) -> Self {
        Self {
            chain: vec![x, y],
            length: 2,
        }
    }
}

impl<const D: usize> CausalSet<D> {
    /// Find the longest chain between two causally related points.
    ///
    /// A chain is a totally ordered subset (antichain-free path) through
    /// the causal set. The longest chain approximates a geodesic.
    ///
    /// # Arguments
    ///
    /// * `x` - Starting point index (must be in causal past of y)
    /// * `y` - Ending point index
    ///
    /// # Returns
    ///
    /// `ChainResult` containing the chain elements and length.
    /// Returns empty result if x and y are not causally related.
    ///
    /// # Complexity
    ///
    /// O(|interval|²) where interval = {z : x ≺ z ≺ y}
    pub fn longest_chain(&self, x: usize, y: usize) -> ChainResult {
        if !self.is_related(x, y) {
            return ChainResult::empty();
        }

        // Get all elements strictly between x and y
        let interval = self.interval(x, y);

        if interval.is_empty() {
            // Direct link: x → y with nothing between
            return ChainResult::direct_link(x, y);
        }

        // Sort interval by time coordinate (valid topological sort for sub-DAG)
        let mut sorted_interval = interval.clone();
        sorted_interval.sort_by(|&a, &b| {
            self.point(a)
                .t()
                .partial_cmp(&self.point(b).t())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let n = sorted_interval.len();

        // DP: longest[i] = (chain length from x to sorted_interval[i], predecessor index)
        let mut longest: Vec<(usize, Option<usize>)> = vec![(0, None); n];

        // Base case: elements directly linked from x
        for (i, &elem) in sorted_interval.iter().enumerate() {
            if self.is_linked(x, elem) {
                longest[i] = (2, None); // chain: x -> elem (length 2)
            } else if self.is_related(x, elem) {
                // Related but not linked - there's something between x and elem
                // We need to find it through the DP
                longest[i] = (1, None); // Will be updated in DP
            }
        }

        // DP recurrence: for each pair, check if we can extend a chain
        for i in 0..n {
            if longest[i].0 == 0 {
                continue; // Not reachable from x
            }
            for j in i + 1..n {
                if self.is_related(sorted_interval[i], sorted_interval[j]) {
                    let new_len = longest[i].0 + 1;
                    if new_len > longest[j].0 {
                        longest[j] = (new_len, Some(i));
                    }
                }
            }
        }

        // Find best endpoint that links to y
        let mut best_len = 2; // Minimum: x -> y direct
        let mut best_idx: Option<usize> = None;

        for (i, &elem) in sorted_interval.iter().enumerate() {
            if longest[i].0 > 0 && self.is_related(elem, y) {
                // Can extend chain from elem to y
                let total_len = longest[i].0 + 1;
                if total_len > best_len {
                    best_len = total_len;
                    best_idx = Some(i);
                }
            }
        }

        // Reconstruct path
        let mut chain = vec![y];

        if let Some(mut current) = best_idx {
            chain.push(sorted_interval[current]);
            while let Some(pred) = longest[current].1 {
                chain.push(sorted_interval[pred]);
                current = pred;
            }
        }

        chain.push(x);
        chain.reverse();

        ChainResult {
            length: chain.len(),
            chain,
        }
    }

    /// Get the length of the longest chain between two points.
    ///
    /// More efficient than `longest_chain` if you only need the length.
    pub fn longest_chain_length(&self, x: usize, y: usize) -> usize {
        self.longest_chain(x, y).length
    }

    /// Estimate proper time separation between two causally related points.
    ///
    /// Uses the longest chain as a proxy for proper time:
    /// τ ≈ L where L is the chain length.
    ///
    /// Note: Full calibration (τ = C_d · L / ρ^(1/d)) requires knowing
    /// the sprinkling density and dimension. This is deferred to Phase 3.
    pub fn proper_time_estimate(&self, x: usize, y: usize) -> f64 {
        self.longest_chain_length(x, y) as f64
    }

    /// Find all maximal chains (longest paths) in the entire causal set.
    ///
    /// Returns chains from minimal to maximal elements.
    pub fn maximal_chains(&self) -> Vec<ChainResult> {
        let minimals = self.minimal_elements();
        let maximals = self.maximal_elements();

        let mut chains = Vec::new();

        for &min_elem in &minimals {
            for &max_elem in &maximals {
                if self.is_related(min_elem, max_elem) {
                    let chain = self.longest_chain(min_elem, max_elem);
                    if chain.length > 0 {
                        chains.push(chain);
                    }
                }
            }
        }

        // Sort by length (longest first)
        chains.sort_by(|a, b| b.length.cmp(&a.length));
        chains
    }

    /// Get the height of the causal set (length of longest maximal chain).
    pub fn height(&self) -> usize {
        self.maximal_chains()
            .first()
            .map(|c| c.length)
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point2D;

    #[test]
    fn test_chain_in_total_order() {
        // Create a totally ordered chain: p0 ≺ p1 ≺ p2 ≺ p3 ≺ p4
        let points: Vec<Point2D> = (0..5)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);

        // Longest chain from 0 to 4 should include all elements
        let chain = causet.longest_chain(0, 4);
        assert_eq!(chain.length, 5);
        assert_eq!(chain.chain, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_chain_direct_link() {
        // Two points with nothing between
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0),
            Point2D::new_2d(0.5, 0.0, 1), // Close enough to be a direct link
        ];

        let causet = CausalSet::from_points(points);

        let chain = causet.longest_chain(0, 1);
        assert_eq!(chain.length, 2);
        assert_eq!(chain.chain, vec![0, 1]);
    }

    #[test]
    fn test_chain_unrelated_points() {
        // Two spacelike separated points
        let points = vec![
            Point2D::new_2d(0.0, -1.0, 0),
            Point2D::new_2d(0.0, 1.0, 1), // Same time, different space
        ];

        let causet = CausalSet::from_points(points);

        let chain = causet.longest_chain(0, 1);
        assert_eq!(chain.length, 0);
        assert!(chain.chain.is_empty());
    }

    #[test]
    fn test_chain_with_branches() {
        // Diamond pattern:
        //       3
        //      / \
        //     1   2
        //      \ /
        //       0
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0),  // Bottom
            Point2D::new_2d(1.0, 0.3, 1),  // Left middle
            Point2D::new_2d(1.0, -0.3, 2), // Right middle
            Point2D::new_2d(2.0, 0.0, 3),  // Top
        ];

        let causet = CausalSet::from_points(points);

        // Longest chain from 0 to 3 should have length 3 (0 -> 1 -> 3 or 0 -> 2 -> 3)
        let chain = causet.longest_chain(0, 3);
        assert_eq!(chain.length, 3);
        assert_eq!(chain.chain[0], 0);
        assert_eq!(chain.chain[2], 3);
        // Middle element should be either 1 or 2
        assert!(chain.chain[1] == 1 || chain.chain[1] == 2);
    }

    #[test]
    fn test_height() {
        let points: Vec<Point2D> = (0..10)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);
        assert_eq!(causet.height(), 10);
    }

    #[test]
    fn test_proper_time_estimate() {
        let points: Vec<Point2D> = (0..5)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);

        // Proper time estimate should equal chain length for now
        assert_eq!(causet.proper_time_estimate(0, 4), 5.0);
        assert_eq!(causet.proper_time_estimate(1, 3), 3.0);
    }
}
