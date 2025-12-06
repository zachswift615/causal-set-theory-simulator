//! # causet-core
//!
//! Core data structures for causal set quantum gravity simulations.
//!
//! This crate provides the fundamental building blocks:
//! - [`SpacetimePoint`] — Events in Lorentzian manifolds
//! - [`CausalDiamond`] — Alexandrov intervals for sprinkling
//! - [`Sprinkler`] — Lorentz-invariant Poisson point process
//! - [`CausalSet`] — The discrete spacetime structure
//!
//! # Quick Start
//!
//! ```rust
//! use causet_core::prelude::*;
//!
//! // Create a causal diamond in 2D Minkowski spacetime
//! let diamond = CausalDiamond::<2>::symmetric(4.0);
//!
//! // Sprinkle points with reproducible seed
//! let mut sprinkler = Sprinkler::with_seed(42, 50.0);
//! let result = sprinkler.sprinkle_diamond(&diamond);
//!
//! // Build the causal set
//! let causet = CausalSet::from_sprinkling(result);
//!
//! // Measure the ordering fraction (should be ~0.5 for 2D)
//! println!("Ordering fraction: {:.4}", causet.ordering_fraction());
//! println!("Elements: {}", causet.len());
//! println!("Relations: {}", causet.causal_matrix().relation_count());
//! println!("Links: {}", causet.link_matrix().link_count());
//! ```
//!
//! # Theoretical Background
//!
//! Causal set theory proposes that spacetime is fundamentally discrete,
//! with geometry emerging from causal structure. The key insight:
//!
//! > **Causal structure + cardinality ≈ spacetime geometry**
//!
//! The ordering fraction f = 2R/[N(N-1)] relates directly to spacetime dimension
//! via the Myrheim-Meyer formula. For Minkowski spacetime:
//! - f(2D) = 0.500
//! - f(3D) ≈ 0.424
//! - f(4D) = 0.333
//!
//! # References
//!
//! - Bombelli, Lee, Meyer, Sorkin (1987) "Space-time as a causal set"
//! - Surya (2019) "The causal set approach to quantum gravity"

pub mod causal_set;
pub mod point;
pub mod relations;
pub mod sprinkling;

/// Convenient re-exports for common usage.
pub mod prelude {
    pub use crate::causal_set::{CausalSet, CausalSetBuilder, CausalSetMetadata};
    pub use crate::point::{Point2D, Point3D, Point4D, SpacetimePoint};
    pub use crate::relations::{
        CausalMatrix, CausalMatrixBuilder, CausalStatistics, LinkMatrix,
    };
    pub use crate::sprinkling::{CausalDiamond, Sprinkler, SprinklingConfig, SprinklingResult};
}

// Re-export main types at crate root
pub use causal_set::CausalSet;
pub use point::SpacetimePoint;
pub use relations::CausalMatrix;
pub use sprinkling::{CausalDiamond, Sprinkler};
