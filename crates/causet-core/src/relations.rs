//! Causal relation computation and storage.
//!
//! This module handles the O(N²) computation of which pairs of points are
//! causally related, along with efficient storage formats.
//!
//! # Key Concepts
//! - **Causal matrix**: C[i,j] = 1 iff point i ≺ point j (i is in causal past of j)
//! - **Link matrix**: L[i,j] = 1 iff i ≺ j and no k exists with i ≺ k ≺ j
//! - Links are the "nearest neighbor" causal relations (transitive reduction)

use crate::point::SpacetimePoint;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use tracing::{debug, info, instrument};

/// Compressed storage for the causal relation matrix.
///
/// Since the matrix is strictly upper triangular (after time-sorting),
/// we only store the upper triangle as a bit vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CausalMatrix {
    /// Number of elements in the causal set
    n: usize,

    /// Bit-packed upper triangular matrix
    /// Index (i,j) with i < j maps to bit position: i*n - i*(i+1)/2 + (j-i-1)
    bits: Vec<u64>,

    /// Number of causal relations (pairs where i ≺ j)
    relation_count: usize,
}

impl CausalMatrix {
    /// Create a new empty causal matrix for n elements.
    pub fn new(n: usize) -> Self {
        // Number of pairs: n*(n-1)/2
        let n_pairs = n * (n.saturating_sub(1)) / 2;
        let n_words = (n_pairs + 63) / 64;

        Self {
            n,
            bits: vec![0u64; n_words],
            relation_count: 0,
        }
    }

    /// Get the number of elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    /// Check if empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Map (i, j) with i < j to bit index.
    #[inline]
    fn index(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < j && j < self.n);
        // Row i has entries for columns i+1, i+2, ..., n-1
        // Number of entries before row i: sum_{k=0}^{i-1} (n-1-k) = i*n - i*(i+1)/2
        // Position within row i: j - i - 1
        i * self.n - i * (i + 1) / 2 + (j - i - 1)
    }

    /// Set the relation i ≺ j (i is in causal past of j).
    ///
    /// Requires i < j (points must be time-sorted).
    #[inline]
    pub fn set_relation(&mut self, i: usize, j: usize) {
        debug_assert!(i < j);
        let idx = self.index(i, j);
        let word = idx / 64;
        let bit = idx % 64;
        if self.bits[word] & (1u64 << bit) == 0 {
            self.bits[word] |= 1u64 << bit;
            self.relation_count += 1;
        }
    }

    /// Check if i ≺ j.
    #[inline]
    pub fn is_related(&self, i: usize, j: usize) -> bool {
        if i >= j {
            return false;
        }
        let idx = self.index(i, j);
        let word = idx / 64;
        let bit = idx % 64;
        (self.bits[word] & (1u64 << bit)) != 0
    }

    /// Get the total number of causal relations.
    #[inline]
    pub fn relation_count(&self) -> usize {
        self.relation_count
    }

    /// Compute the ordering fraction f = 2R / [N(N-1)].
    ///
    /// This is the key observable for dimension estimation.
    pub fn ordering_fraction(&self) -> f64 {
        if self.n <= 1 {
            return 0.0;
        }
        2.0 * self.relation_count as f64 / (self.n * (self.n - 1)) as f64
    }

    /// Iterate over all related pairs (i, j) where i ≺ j.
    pub fn related_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.n).flat_map(move |i| {
            (i + 1..self.n).filter_map(move |j| {
                if self.is_related(i, j) {
                    Some((i, j))
                } else {
                    None
                }
            })
        })
    }

    /// Get all elements in the causal past of j (i.e., all i where i ≺ j).
    pub fn past_of(&self, j: usize) -> Vec<usize> {
        (0..j).filter(|&i| self.is_related(i, j)).collect()
    }

    /// Get all elements in the causal future of i (i.e., all j where i ≺ j).
    pub fn future_of(&self, i: usize) -> Vec<usize> {
        (i + 1..self.n)
            .filter(|&j| self.is_related(i, j))
            .collect()
    }
}

/// Builder for constructing causal matrices from point sets.
pub struct CausalMatrixBuilder {
    parallel: bool,
}

impl Default for CausalMatrixBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl CausalMatrixBuilder {
    pub fn new() -> Self {
        Self { parallel: true }
    }

    /// Enable or disable parallel computation.
    pub fn parallel(mut self, parallel: bool) -> Self {
        self.parallel = parallel;
        self
    }

    /// Build the causal matrix from a set of points.
    ///
    /// Points are first sorted by time coordinate, then all pairs are checked.
    /// Complexity: O(N²) comparisons, O(N log N) sort.
    #[instrument(skip(self, points), fields(n = points.len()))]
    pub fn build<const D: usize>(&self, points: &[SpacetimePoint<D>]) -> CausalMatrixResult {
        let n = points.len();
        info!(n = n, "Building causal matrix");

        // Sort points by time coordinate
        let mut sorted_indices: Vec<usize> = (0..n).collect();
        sorted_indices.sort_by(|&a, &b| {
            points[a]
                .t()
                .partial_cmp(&points[b].t())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Create sorted point references
        let sorted_points: Vec<&SpacetimePoint<D>> =
            sorted_indices.iter().map(|&i| &points[i]).collect();

        let mut matrix = CausalMatrix::new(n);

        if self.parallel && n > 100 {
            // Parallel computation for large sets
            let relations: Vec<(usize, usize)> = (0..n)
                .into_par_iter()
                .flat_map(|i| {
                    (i + 1..n)
                        .filter_map(|j| {
                            if sorted_points[i].precedes(sorted_points[j]) {
                                Some((i, j))
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect();

            for (i, j) in relations {
                matrix.set_relation(i, j);
            }
        } else {
            // Sequential computation for small sets
            for i in 0..n {
                for j in i + 1..n {
                    if sorted_points[i].precedes(sorted_points[j]) {
                        matrix.set_relation(i, j);
                    }
                }
            }
        }

        debug!(
            relations = matrix.relation_count(),
            ordering_fraction = matrix.ordering_fraction(),
            "Causal matrix complete"
        );

        CausalMatrixResult {
            matrix,
            sorted_indices,
        }
    }
}

/// Result of building a causal matrix.
#[derive(Debug)]
pub struct CausalMatrixResult {
    /// The computed causal matrix (indices refer to sorted order)
    pub matrix: CausalMatrix,

    /// Mapping from sorted index to original point index
    pub sorted_indices: Vec<usize>,
}

impl CausalMatrixResult {
    /// Map a sorted index back to the original point index.
    pub fn original_index(&self, sorted_idx: usize) -> usize {
        self.sorted_indices[sorted_idx]
    }
}

/// Compute the transitive reduction (link matrix) from a causal matrix.
///
/// A link exists between i and j iff:
/// - i ≺ j (they are causally related)
/// - There is no k such that i ≺ k ≺ j (the relation is "direct")
///
/// Complexity: O(N³) in worst case, but typically much better for manifold-like causets.
#[instrument(skip(matrix), fields(n = matrix.len()))]
pub fn compute_links(matrix: &CausalMatrix) -> LinkMatrix {
    let n = matrix.len();
    info!(n = n, "Computing transitive reduction (links)");

    let mut links = LinkMatrix::new(n);
    let mut link_count = 0usize;

    // For each related pair, check if there's an intermediate element
    for i in 0..n {
        for j in i + 1..n {
            if !matrix.is_related(i, j) {
                continue;
            }

            // Check for any k with i ≺ k ≺ j
            let mut is_link = true;
            for k in i + 1..j {
                if matrix.is_related(i, k) && matrix.is_related(k, j) {
                    is_link = false;
                    break;
                }
            }

            if is_link {
                links.set_link(i, j);
                link_count += 1;
            }
        }
    }

    debug!(link_count = link_count, "Transitive reduction complete");
    links
}

/// Storage for link relations (transitive reduction of causal matrix).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinkMatrix {
    n: usize,
    bits: Vec<u64>,
    link_count: usize,
}

impl LinkMatrix {
    pub fn new(n: usize) -> Self {
        let n_pairs = n * (n.saturating_sub(1)) / 2;
        let n_words = (n_pairs + 63) / 64;

        Self {
            n,
            bits: vec![0u64; n_words],
            link_count: 0,
        }
    }

    #[inline]
    fn index(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < j && j < self.n);
        i * self.n - i * (i + 1) / 2 + (j - i - 1)
    }

    #[inline]
    pub fn set_link(&mut self, i: usize, j: usize) {
        debug_assert!(i < j);
        let idx = self.index(i, j);
        let word = idx / 64;
        let bit = idx % 64;
        if self.bits[word] & (1u64 << bit) == 0 {
            self.bits[word] |= 1u64 << bit;
            self.link_count += 1;
        }
    }

    #[inline]
    pub fn is_linked(&self, i: usize, j: usize) -> bool {
        if i >= j {
            return false;
        }
        let idx = self.index(i, j);
        let word = idx / 64;
        let bit = idx % 64;
        (self.bits[word] & (1u64 << bit)) != 0
    }

    pub fn link_count(&self) -> usize {
        self.link_count
    }

    /// Get direct children (future links) of element i.
    pub fn children(&self, i: usize) -> Vec<usize> {
        (i + 1..self.n)
            .filter(|&j| self.is_linked(i, j))
            .collect()
    }

    /// Get direct parents (past links) of element j.
    pub fn parents(&self, j: usize) -> Vec<usize> {
        (0..j).filter(|&i| self.is_linked(i, j)).collect()
    }

    /// Iterate over all links as (parent, child) pairs.
    pub fn links(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.n).flat_map(move |i| {
            (i + 1..self.n).filter_map(move |j| {
                if self.is_linked(i, j) {
                    Some((i, j))
                } else {
                    None
                }
            })
        })
    }
}

/// Statistics about the causal structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CausalStatistics {
    pub n_elements: usize,
    pub n_relations: usize,
    pub n_links: usize,
    pub ordering_fraction: f64,
    pub link_density: f64,
}

impl CausalStatistics {
    pub fn compute(matrix: &CausalMatrix, links: &LinkMatrix) -> Self {
        let n = matrix.len();
        let n_pairs = if n > 1 { n * (n - 1) / 2 } else { 1 };

        Self {
            n_elements: n,
            n_relations: matrix.relation_count(),
            n_links: links.link_count(),
            ordering_fraction: matrix.ordering_fraction(),
            link_density: links.link_count() as f64 / n_pairs as f64,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point2D;

    #[test]
    fn test_causal_matrix_basic() {
        let mut matrix = CausalMatrix::new(5);

        matrix.set_relation(0, 2);
        matrix.set_relation(0, 4);
        matrix.set_relation(1, 3);

        assert!(matrix.is_related(0, 2));
        assert!(matrix.is_related(0, 4));
        assert!(matrix.is_related(1, 3));
        assert!(!matrix.is_related(0, 1));
        assert!(!matrix.is_related(2, 4));

        assert_eq!(matrix.relation_count(), 3);
    }

    #[test]
    fn test_build_causal_matrix() {
        // Create a simple chain: p0 ≺ p1 ≺ p2
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0),
            Point2D::new_2d(1.0, 0.0, 1),
            Point2D::new_2d(2.0, 0.0, 2),
        ];

        let result = CausalMatrixBuilder::new().parallel(false).build(&points);

        // All points should be causally related in a chain
        assert!(result.matrix.is_related(0, 1));
        assert!(result.matrix.is_related(1, 2));
        assert!(result.matrix.is_related(0, 2)); // Transitive

        assert_eq!(result.matrix.relation_count(), 3);
    }

    #[test]
    fn test_transitive_reduction() {
        // Create a chain: p0 ≺ p1 ≺ p2
        // Links should be: p0 → p1 → p2 (not p0 → p2)
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0),
            Point2D::new_2d(1.0, 0.0, 1),
            Point2D::new_2d(2.0, 0.0, 2),
        ];

        let result = CausalMatrixBuilder::new().parallel(false).build(&points);
        let links = compute_links(&result.matrix);

        assert!(links.is_linked(0, 1));
        assert!(links.is_linked(1, 2));
        assert!(!links.is_linked(0, 2)); // Transitive relation, not a link

        assert_eq!(links.link_count(), 2);
    }

    #[test]
    fn test_spacelike_separated() {
        // Two points at same time (spacelike separated)
        let points = vec![
            Point2D::new_2d(0.0, -1.0, 0),
            Point2D::new_2d(0.0, 1.0, 1),
        ];

        let result = CausalMatrixBuilder::new().parallel(false).build(&points);

        assert_eq!(result.matrix.relation_count(), 0);
    }

    #[test]
    fn test_ordering_fraction() {
        // For 2D Minkowski spacetime, ordering fraction should be ~0.5
        // Create points that should give roughly this fraction

        // Simple diamond-like pattern
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0), // Bottom
            Point2D::new_2d(1.0, 0.5, 1), // Middle-right (in future of 0)
            Point2D::new_2d(1.0, -0.5, 2), // Middle-left (in future of 0)
            Point2D::new_2d(2.0, 0.0, 3), // Top (in future of all)
        ];

        let result = CausalMatrixBuilder::new().parallel(false).build(&points);

        // Relations: 0≺1, 0≺2, 0≺3, 1≺3, 2≺3 = 5 relations
        // 1 and 2 are spacelike (Δt=0)
        // Total pairs: 4*3/2 = 6
        // f = 2*5/12 ≈ 0.833 (this is a highly connected set)

        assert!(result.matrix.ordering_fraction() > 0.0);
        assert!(result.matrix.ordering_fraction() <= 1.0);
    }
}
