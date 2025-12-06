//! Spacetime point representations for causal set theory.
//!
//! This module defines the fundamental `SpacetimePoint` type representing
//! events in Lorentzian manifolds. Points carry coordinates and metadata
//! for tracking their origin and properties.

use nalgebra::{SVector, Vector2, Vector3, Vector4};
use serde::{Deserialize, Serialize};

/// A point in d-dimensional Minkowski spacetime.
///
/// Coordinates use the (-,+,+,...) signature convention where:
/// - `coords[0]` is the time coordinate t
/// - `coords[1..]` are spatial coordinates x, y, z, ...
///
/// The Minkowski metric is ds² = -dt² + dx² + dy² + dz² + ...
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpacetimePoint<const D: usize> {
    /// Coordinates in the embedding manifold [t, x, y, z, ...]
    pub coords: SVector<f64, D>,

    /// Unique identifier assigned during sprinkling
    pub id: usize,

    /// Index of the sprinkling batch that created this point (for reproducibility)
    pub batch_id: u64,
}

impl<const D: usize> SpacetimePoint<D> {
    /// Create a new spacetime point with given coordinates.
    pub fn new(coords: SVector<f64, D>, id: usize) -> Self {
        Self {
            coords,
            id,
            batch_id: 0,
        }
    }

    /// Create a new spacetime point with batch tracking.
    pub fn with_batch(coords: SVector<f64, D>, id: usize, batch_id: u64) -> Self {
        Self {
            coords,
            id,
            batch_id,
        }
    }

    /// Get the time coordinate t.
    #[inline]
    pub fn t(&self) -> f64 {
        self.coords[0]
    }

    /// Get the spatial coordinates as a slice.
    #[inline]
    pub fn spatial(&self) -> &[f64] {
        &self.coords.as_slice()[1..]
    }

    /// Compute the Minkowski interval squared between two points.
    ///
    /// Returns Δs² = -(Δt)² + |Δx|²
    /// - Negative (timelike): causal connection possible
    /// - Zero (null/lightlike): on light cone
    /// - Positive (spacelike): no causal connection
    #[inline]
    pub fn interval_squared(&self, other: &Self) -> f64 {
        let dt = other.coords[0] - self.coords[0];
        let dx_squared: f64 = (1..D)
            .map(|i| {
                let dx = other.coords[i] - self.coords[i];
                dx * dx
            })
            .sum();

        -dt * dt + dx_squared
    }

    /// Check if this point is in the causal past of another point.
    ///
    /// Returns true iff:
    /// 1. The interval is strictly timelike (Δs² < 0)
    /// 2. The other point is in the future (t_other > t_self)
    ///
    /// Note: We use STRICT inequality (>) to exclude null-separated points.
    /// This matches the Myrheim-Meyer formula convention for ordering fraction.
    #[inline]
    pub fn precedes(&self, other: &Self) -> bool {
        let dt = other.coords[0] - self.coords[0];
        if dt <= 0.0 {
            return false;
        }

        let dx_squared: f64 = (1..D)
            .map(|i| {
                let dx = other.coords[i] - self.coords[i];
                dx * dx
            })
            .sum();

        // Strictly timelike: dt² > |dx|² (excludes null/lightlike separation)
        dt * dt > dx_squared
    }
}

/// Type alias for 2D Minkowski spacetime (1+1 dimensions).
pub type Point2D = SpacetimePoint<2>;

/// Type alias for 3D Minkowski spacetime (2+1 dimensions).
pub type Point3D = SpacetimePoint<3>;

/// Type alias for 4D Minkowski spacetime (3+1 dimensions).
pub type Point4D = SpacetimePoint<4>;

// Convenience constructors for common dimensions
impl SpacetimePoint<2> {
    /// Create a 2D point from (t, x) coordinates.
    pub fn new_2d(t: f64, x: f64, id: usize) -> Self {
        Self::new(Vector2::new(t, x), id)
    }
}

impl SpacetimePoint<3> {
    /// Create a 3D point from (t, x, y) coordinates.
    pub fn new_3d(t: f64, x: f64, y: f64, id: usize) -> Self {
        Self::new(Vector3::new(t, x, y), id)
    }
}

impl SpacetimePoint<4> {
    /// Create a 4D point from (t, x, y, z) coordinates.
    pub fn new_4d(t: f64, x: f64, y: f64, z: f64, id: usize) -> Self {
        Self::new(Vector4::new(t, x, y, z), id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_timelike_interval() {
        // Two points with Δt = 2, Δx = 1 → Δs² = -4 + 1 = -3 (timelike)
        let p1 = Point2D::new_2d(0.0, 0.0, 0);
        let p2 = Point2D::new_2d(2.0, 1.0, 1);

        assert_relative_eq!(p1.interval_squared(&p2), -3.0);
        assert!(p1.precedes(&p2));
        assert!(!p2.precedes(&p1));
    }

    #[test]
    fn test_spacelike_interval() {
        // Two points with Δt = 1, Δx = 2 → Δs² = -1 + 4 = 3 (spacelike)
        let p1 = Point2D::new_2d(0.0, 0.0, 0);
        let p2 = Point2D::new_2d(1.0, 2.0, 1);

        assert_relative_eq!(p1.interval_squared(&p2), 3.0);
        assert!(!p1.precedes(&p2));
        assert!(!p2.precedes(&p1));
    }

    #[test]
    fn test_null_interval() {
        // Two points on the light cone: Δt = 1, Δx = 1 → Δs² = 0
        let p1 = Point2D::new_2d(0.0, 0.0, 0);
        let p2 = Point2D::new_2d(1.0, 1.0, 1);

        assert_relative_eq!(p1.interval_squared(&p2), 0.0, epsilon = 1e-10);
        // Null-separated points are NOT causally related (strict timelike required)
        assert!(!p1.precedes(&p2));
    }

    #[test]
    fn test_4d_interval() {
        // 4D: Δt = 3, Δx = 1, Δy = 1, Δz = 1 → Δs² = -9 + 3 = -6 (timelike)
        let p1 = Point4D::new_4d(0.0, 0.0, 0.0, 0.0, 0);
        let p2 = Point4D::new_4d(3.0, 1.0, 1.0, 1.0, 1);

        assert_relative_eq!(p1.interval_squared(&p2), -6.0);
        assert!(p1.precedes(&p2));
    }

    #[test]
    fn test_same_time_not_causal() {
        // Points at same time cannot be causally related
        let p1 = Point2D::new_2d(1.0, 0.0, 0);
        let p2 = Point2D::new_2d(1.0, 5.0, 1);

        assert!(!p1.precedes(&p2));
        assert!(!p2.precedes(&p1));
    }
}
