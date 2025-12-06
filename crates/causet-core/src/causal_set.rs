//! The main CausalSet data structure.
//!
//! A causal set is a locally finite partially ordered set (poset) where:
//! 1. The order relation ≺ is transitive and irreflexive
//! 2. Every interval [x,y] = {z : x ≺ z ≺ y} is finite
//!
//! The fundamental conjecture: causal structure + cardinality ≈ spacetime geometry.

use crate::point::SpacetimePoint;
use crate::relations::{
    compute_links, CausalMatrix, CausalMatrixBuilder, CausalStatistics, LinkMatrix,
};
use crate::spacetime::Spacetime;
use crate::sprinkling::{GenericSprinklingResult, SprinklingResult};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use serde::{Deserialize, Serialize};
use tracing::{debug, info, instrument};

/// A causal set: discrete spacetime structure.
///
/// Stores both the full causal matrix (for fast queries) and the DAG of links
/// (for graph algorithms and visualization).
#[derive(Debug)]
pub struct CausalSet<const D: usize> {
    /// The spacetime points (sorted by time coordinate)
    points: Vec<SpacetimePoint<D>>,

    /// Mapping from sorted index to original sprinkling index
    original_indices: Vec<usize>,

    /// Full causal relation matrix
    causal_matrix: CausalMatrix,

    /// Link matrix (transitive reduction)
    link_matrix: LinkMatrix,

    /// DAG representation using petgraph (edges are links)
    dag: DiGraph<usize, ()>,

    /// Node indices in the DAG (maps sorted point index to NodeIndex)
    node_indices: Vec<NodeIndex>,

    /// Metadata for reproducibility
    metadata: CausalSetMetadata,
}

/// Metadata for tracking causal set provenance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CausalSetMetadata {
    pub dimensions: usize,
    pub sprinkling_seed: Option<u64>,
    pub sprinkling_density: Option<f64>,
    pub diamond_proper_time: Option<f64>,
    pub diamond_volume: Option<f64>,
}

impl<const D: usize> CausalSet<D> {
    /// Build a causal set from a sprinkling result.
    #[instrument(skip(result), fields(n = result.points.len()))]
    pub fn from_sprinkling(result: SprinklingResult<D>) -> Self {
        info!(
            n = result.points.len(),
            density = result.config.density,
            "Building causal set from sprinkling"
        );

        let metadata = CausalSetMetadata {
            dimensions: D,
            sprinkling_seed: Some(result.config.seed),
            sprinkling_density: Some(result.config.density),
            diamond_proper_time: Some(result.diamond.proper_time),
            diamond_volume: Some(result.diamond.volume()),
        };

        Self::from_points_with_metadata(result.points, metadata)
    }

    /// Build a causal set from a vector of points.
    pub fn from_points(points: Vec<SpacetimePoint<D>>) -> Self {
        let metadata = CausalSetMetadata {
            dimensions: D,
            sprinkling_seed: None,
            sprinkling_density: None,
            diamond_proper_time: None,
            diamond_volume: None,
        };

        Self::from_points_with_metadata(points, metadata)
    }

    /// Build a causal set from a generic sprinkling result using a specific spacetime.
    ///
    /// This is the Phase 4 method for curved spacetime sprinkling. It uses the
    /// spacetime's `causally_precedes` method to determine causal relations,
    /// which may differ from flat Minkowski spacetime.
    ///
    /// # Arguments
    ///
    /// * `result` - The sprinkling result containing points
    /// * `spacetime` - The spacetime geometry for determining causal relations
    ///
    /// # Example
    ///
    /// ```rust
    /// use causet_core::prelude::*;
    ///
    /// // Sprinkle into de Sitter
    /// let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
    /// let result = sprinkle_spacetime(&ds, 500, 42);
    ///
    /// // Build causal set using de Sitter causal structure
    /// let causet = CausalSet::from_generic_sprinkling(result, &ds);
    /// ```
    #[instrument(skip(result, spacetime), fields(n = result.points.len(), spacetime = %result.spacetime_name))]
    pub fn from_generic_sprinkling<S: Spacetime<D>>(
        result: GenericSprinklingResult<D>,
        spacetime: &S,
    ) -> Self {
        debug_assert_eq!(
            result.spacetime_name,
            spacetime.name(),
            "Spacetime mismatch: sprinkled in '{}' but building with '{}'",
            result.spacetime_name,
            spacetime.name()
        );

        info!(
            n = result.points.len(),
            spacetime = %result.spacetime_name,
            ricci_scalar = result.ricci_scalar,
            "Building causal set from generic sprinkling"
        );

        let metadata = CausalSetMetadata {
            dimensions: D,
            sprinkling_seed: Some(result.seed),
            sprinkling_density: Some(result.points.len() as f64 / result.volume),
            diamond_proper_time: None,
            diamond_volume: Some(result.volume),
        };

        Self::from_points_with_spacetime(result.points, spacetime, metadata)
    }

    /// Build a causal set from points using a specific spacetime for causal relations.
    fn from_points_with_spacetime<S: Spacetime<D>>(
        mut points: Vec<SpacetimePoint<D>>,
        spacetime: &S,
        metadata: CausalSetMetadata,
    ) -> Self {
        let n = points.len();

        // Sort points by time coordinate
        let mut sorted_indices: Vec<usize> = (0..n).collect();
        sorted_indices.sort_by(|&a, &b| {
            points[a]
                .t()
                .partial_cmp(&points[b].t())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Reorder points to match sorted order
        let sorted_points: Vec<SpacetimePoint<D>> = sorted_indices
            .iter()
            .map(|&i| points[i].clone())
            .collect();
        points = sorted_points;
        let original_indices = sorted_indices;

        // Build causal matrix using spacetime's causal relation
        let mut causal_matrix = CausalMatrix::new(n);
        for i in 0..n {
            for j in i + 1..n {
                if spacetime.causally_precedes(&points[i].coords, &points[j].coords) {
                    causal_matrix.set_relation(i, j);
                }
            }
        }

        // Compute links
        let link_matrix = compute_links(&causal_matrix);

        // Build petgraph DAG from links
        let mut dag = DiGraph::new();
        let node_indices: Vec<NodeIndex> = (0..n).map(|i| dag.add_node(i)).collect();

        for (i, j) in link_matrix.links() {
            dag.add_edge(node_indices[i], node_indices[j], ());
        }

        debug!(
            n = n,
            relations = causal_matrix.relation_count(),
            links = link_matrix.link_count(),
            "Causal set construction complete (with spacetime)"
        );

        Self {
            points,
            original_indices,
            causal_matrix,
            link_matrix,
            dag,
            node_indices,
            metadata,
        }
    }

    fn from_points_with_metadata(
        mut points: Vec<SpacetimePoint<D>>,
        metadata: CausalSetMetadata,
    ) -> Self {
        let n = points.len();

        // Build causal matrix (sorts points internally)
        let matrix_result = CausalMatrixBuilder::new().build(&points);
        let causal_matrix = matrix_result.matrix;
        let original_indices = matrix_result.sorted_indices.clone();

        // Reorder points to match sorted order
        let sorted_points: Vec<SpacetimePoint<D>> = matrix_result
            .sorted_indices
            .iter()
            .map(|&i| points[i].clone())
            .collect();
        points = sorted_points;

        // Compute links
        let link_matrix = compute_links(&causal_matrix);

        // Build petgraph DAG from links
        let mut dag = DiGraph::new();
        let node_indices: Vec<NodeIndex> = (0..n).map(|i| dag.add_node(i)).collect();

        for (i, j) in link_matrix.links() {
            dag.add_edge(node_indices[i], node_indices[j], ());
        }

        debug!(
            n = n,
            relations = causal_matrix.relation_count(),
            links = link_matrix.link_count(),
            "Causal set construction complete"
        );

        Self {
            points,
            original_indices,
            causal_matrix,
            link_matrix,
            dag,
            node_indices,
            metadata,
        }
    }

    /// Get the number of elements.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Check if the causal set is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Get a point by its sorted index.
    pub fn point(&self, idx: usize) -> &SpacetimePoint<D> {
        &self.points[idx]
    }

    /// Get all points (in time-sorted order).
    pub fn points(&self) -> &[SpacetimePoint<D>] {
        &self.points
    }

    /// Check if two elements are causally related (i ≺ j).
    pub fn is_related(&self, i: usize, j: usize) -> bool {
        self.causal_matrix.is_related(i, j)
    }

    /// Check if two elements are linked (direct causal relation).
    pub fn is_linked(&self, i: usize, j: usize) -> bool {
        self.link_matrix.is_linked(i, j)
    }

    /// Get the causal matrix.
    pub fn causal_matrix(&self) -> &CausalMatrix {
        &self.causal_matrix
    }

    /// Get the link matrix.
    pub fn link_matrix(&self) -> &LinkMatrix {
        &self.link_matrix
    }

    /// Get the underlying DAG.
    pub fn dag(&self) -> &DiGraph<usize, ()> {
        &self.dag
    }

    /// Get the ordering fraction (key observable for dimension estimation).
    ///
    /// f = 2R / [N(N-1)] where R is the number of related pairs.
    pub fn ordering_fraction(&self) -> f64 {
        self.causal_matrix.ordering_fraction()
    }

    /// Get the causal past of an element (all i where i ≺ j).
    pub fn past(&self, j: usize) -> Vec<usize> {
        self.causal_matrix.past_of(j)
    }

    /// Get the causal future of an element (all k where j ≺ k).
    pub fn future(&self, j: usize) -> Vec<usize> {
        self.causal_matrix.future_of(j)
    }

    /// Get direct past (parents via links).
    pub fn parents(&self, j: usize) -> Vec<usize> {
        self.dag
            .edges_directed(self.node_indices[j], Direction::Incoming)
            .map(|e| self.dag[e.source()])
            .collect()
    }

    /// Get direct future (children via links).
    pub fn children(&self, i: usize) -> Vec<usize> {
        self.dag
            .edges_directed(self.node_indices[i], Direction::Outgoing)
            .map(|e| self.dag[e.target()])
            .collect()
    }

    /// Get minimal elements (no past).
    pub fn minimal_elements(&self) -> Vec<usize> {
        (0..self.len())
            .filter(|&i| self.past(i).is_empty())
            .collect()
    }

    /// Get maximal elements (no future).
    pub fn maximal_elements(&self) -> Vec<usize> {
        (0..self.len())
            .filter(|&i| self.future(i).is_empty())
            .collect()
    }

    /// Compute the causal interval [x, y] = {z : x ≺ z ≺ y}.
    ///
    /// This is fundamental for the Benincasa-Dowker action.
    pub fn interval(&self, x: usize, y: usize) -> Vec<usize> {
        if !self.is_related(x, y) {
            return vec![];
        }

        (x + 1..y)
            .filter(|&z| self.is_related(x, z) && self.is_related(z, y))
            .collect()
    }

    /// Count the size of the interval [x, y].
    pub fn interval_size(&self, x: usize, y: usize) -> usize {
        self.interval(x, y).len()
    }

    /// Compute statistics about the causal structure.
    pub fn statistics(&self) -> CausalStatistics {
        CausalStatistics::compute(&self.causal_matrix, &self.link_matrix)
    }

    /// Get metadata for reproducibility documentation.
    pub fn metadata(&self) -> &CausalSetMetadata {
        &self.metadata
    }

    /// Get the mapping from sorted index to original sprinkling index.
    /// Useful for correlating with original point data.
    pub fn original_indices(&self) -> &[usize] {
        &self.original_indices
    }

    /// Export full state for research documentation.
    pub fn to_research_json(&self) -> String {
        let stats = self.statistics();

        serde_json::to_string_pretty(&ResearchExport {
            metadata: self.metadata.clone(),
            statistics: stats,
            point_count: self.len(),
            ordering_fraction: self.ordering_fraction(),
        })
        .unwrap_or_default()
    }
}

#[derive(Serialize)]
struct ResearchExport {
    metadata: CausalSetMetadata,
    statistics: CausalStatistics,
    point_count: usize,
    ordering_fraction: f64,
}

/// Builder for creating causal sets with custom options.
pub struct CausalSetBuilder<const D: usize> {
    points: Vec<SpacetimePoint<D>>,
    metadata: CausalSetMetadata,
}

impl<const D: usize> CausalSetBuilder<D> {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            metadata: CausalSetMetadata {
                dimensions: D,
                sprinkling_seed: None,
                sprinkling_density: None,
                diamond_proper_time: None,
                diamond_volume: None,
            },
        }
    }

    pub fn add_point(mut self, point: SpacetimePoint<D>) -> Self {
        self.points.push(point);
        self
    }

    pub fn add_points(mut self, points: impl IntoIterator<Item = SpacetimePoint<D>>) -> Self {
        self.points.extend(points);
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.metadata.sprinkling_seed = Some(seed);
        self
    }

    pub fn build(self) -> CausalSet<D> {
        CausalSet::from_points_with_metadata(self.points, self.metadata)
    }
}

impl<const D: usize> Default for CausalSetBuilder<D> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point::Point2D;
    use crate::sprinkling::{CausalDiamond, Sprinkler};

    #[test]
    fn test_causal_set_from_points() {
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0),
            Point2D::new_2d(1.0, 0.0, 1),
            Point2D::new_2d(2.0, 0.0, 2),
        ];

        let causet = CausalSet::from_points(points);

        assert_eq!(causet.len(), 3);
        assert!(causet.is_related(0, 1));
        assert!(causet.is_related(1, 2));
        assert!(causet.is_related(0, 2));

        // Links should only be direct
        assert!(causet.is_linked(0, 1));
        assert!(causet.is_linked(1, 2));
        assert!(!causet.is_linked(0, 2));
    }

    #[test]
    fn test_interval_computation() {
        // Chain of 5 elements
        let points: Vec<Point2D> = (0..5)
            .map(|i| Point2D::new_2d(i as f64, 0.0, i))
            .collect();

        let causet = CausalSet::from_points(points);

        // Interval [0, 4] should contain {1, 2, 3}
        let interval = causet.interval(0, 4);
        assert_eq!(interval, vec![1, 2, 3]);

        // Interval [0, 2] should contain {1}
        let interval = causet.interval(0, 2);
        assert_eq!(interval, vec![1]);

        // Interval [1, 2] should be empty (adjacent elements)
        let interval = causet.interval(1, 2);
        assert!(interval.is_empty());
    }

    #[test]
    fn test_from_sprinkling() {
        let diamond = CausalDiamond::<2>::symmetric(4.0);
        let mut sprinkler = Sprinkler::with_seed(42, 10.0);
        let result = sprinkler.sprinkle_diamond(&diamond);

        let causet = CausalSet::from_sprinkling(result);

        assert!(causet.len() > 0);
        assert!(causet.ordering_fraction() > 0.0);
        assert!(causet.ordering_fraction() <= 1.0);

        // Metadata should be preserved
        assert_eq!(causet.metadata().sprinkling_seed, Some(42));
        assert_eq!(causet.metadata().sprinkling_density, Some(10.0));
    }

    #[test]
    fn test_dag_structure() {
        let points = vec![
            Point2D::new_2d(0.0, 0.0, 0),
            Point2D::new_2d(1.0, 0.5, 1),
            Point2D::new_2d(1.0, -0.5, 2),
            Point2D::new_2d(2.0, 0.0, 3),
        ];

        let causet = CausalSet::from_points(points);

        // Check DAG edge count matches link count
        assert_eq!(causet.dag().edge_count(), causet.link_matrix().link_count());

        // Element 0 should be minimal
        assert!(causet.minimal_elements().contains(&0));

        // Element 3 should be maximal
        assert!(causet.maximal_elements().contains(&3));
    }
}
