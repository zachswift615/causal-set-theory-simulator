//! Spacetime geometry abstractions for causal set sprinkling.
//!
//! This module provides the [`Spacetime`] trait that defines how to:
//! - Sample points with correct volume weighting
//! - Check causal relations between points
//! - Query geometric properties (curvature, volume)
//!
//! # Implementations
//!
//! - [`Minkowski`] — Flat Minkowski spacetime (causal diamond)
//! - [`DeSitter`] — de Sitter spacetime (conformal coordinates)
//!
//! # Coordinate Convention
//!
//! All spacetime implementations use the convention:
//! - Index 0 is the timelike coordinate
//! - Indices 1..D are spacelike coordinates
//! - Metric signature is (-,+,+,+)

mod de_sitter;
mod minkowski;

pub use de_sitter::DeSitter;
pub use minkowski::Minkowski;

use nalgebra::SVector;
use rand::Rng;

/// Trait defining a spacetime geometry for causal set sprinkling.
///
/// This trait encapsulates all the geometric information needed to:
/// 1. Generate points with the correct volume measure
/// 2. Determine causal relations between points
/// 3. Compute geometric observables
///
/// # Coordinate Convention
///
/// All implementations must follow:
/// - `coords[0]` = timelike coordinate
/// - `coords[1..]` = spacelike coordinates
/// - Metric signature (-,+,+,+)
///
/// # Conformally Flat Spacetimes
///
/// For conformally flat spacetimes (Minkowski, de Sitter, FLRW), the causal
/// structure is identical — only the volume measure differs. Use
/// [`conformal_causal_check`] to share the causal relation logic.
pub trait Spacetime<const D: usize>: Clone + Send + Sync {
    /// Human-readable name for logging and test output.
    ///
    /// Example: "Minkowski-4D", "deSitter-4D (H=0.1)"
    fn name(&self) -> String;

    /// Sample a point uniformly distributed according to the spacetime volume measure.
    ///
    /// The sampling must respect the metric's volume element √(-g).
    /// For curved spacetimes, this typically requires importance sampling.
    fn sample_point(&self, rng: &mut impl Rng) -> SVector<f64, D>;

    /// Check if point p1 causally precedes point p2 (p1 ≺ p2).
    ///
    /// Returns true iff there exists a future-directed causal curve from p1 to p2.
    /// For conformally flat spacetimes, this is equivalent to the Minkowski condition.
    fn causally_precedes(&self, p1: &SVector<f64, D>, p2: &SVector<f64, D>) -> bool;

    /// The Ricci scalar curvature R of the spacetime.
    ///
    /// For maximally symmetric spacetimes:
    /// - Minkowski: R = 0
    /// - de Sitter: R = d(d-1)H² (positive)
    /// - Anti-de Sitter: R = -d(d-1)H² (negative)
    fn ricci_scalar(&self) -> f64;

    /// The proper volume of the sprinkling region.
    ///
    /// This is ∫ √(-g) d^D x over the sampling region.
    fn volume(&self) -> f64;
}

/// Check causal precedence for conformally flat spacetimes.
///
/// For spacetimes with metric ds² = Ω²(x)(-dt² + dx⃗²), the causal structure
/// is identical to Minkowski spacetime. Two points are causally related iff
/// the coordinate separation is timelike in the Minkowski sense.
///
/// This function can be used by any conformally flat spacetime:
/// - Minkowski (Ω = 1)
/// - de Sitter in conformal coordinates (Ω = 1/(H|η|))
/// - FLRW cosmologies (Ω = a(η))
///
/// # Arguments
///
/// * `p1` - First point (potential past point)
/// * `p2` - Second point (potential future point)
///
/// # Returns
///
/// `true` if p1 is in the causal past of p2 (strictly timelike, future-directed)
#[inline]
pub fn conformal_causal_check<const D: usize>(
    p1: &SVector<f64, D>,
    p2: &SVector<f64, D>,
) -> bool {
    let dt = p2[0] - p1[0];
    if dt <= 0.0 {
        return false;
    }

    let dx_sq: f64 = (1..D).map(|i| (p2[i] - p1[i]).powi(2)).sum();

    // Strictly timelike: dt² > |dx|² (excludes null separation)
    dt * dt > dx_sq
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector4;

    #[test]
    fn test_conformal_causal_check_timelike() {
        let p1 = Vector4::new(0.0, 0.0, 0.0, 0.0);
        let p2 = Vector4::new(2.0, 1.0, 0.0, 0.0); // dt=2, |dx|=1 → timelike

        assert!(conformal_causal_check(&p1, &p2));
        assert!(!conformal_causal_check(&p2, &p1)); // Wrong time order
    }

    #[test]
    fn test_conformal_causal_check_spacelike() {
        let p1 = Vector4::new(0.0, 0.0, 0.0, 0.0);
        let p2 = Vector4::new(1.0, 2.0, 0.0, 0.0); // dt=1, |dx|=2 → spacelike

        assert!(!conformal_causal_check(&p1, &p2));
    }

    #[test]
    fn test_conformal_causal_check_null() {
        let p1 = Vector4::new(0.0, 0.0, 0.0, 0.0);
        let p2 = Vector4::new(1.0, 1.0, 0.0, 0.0); // dt=1, |dx|=1 → null

        // Null separation is NOT causal (strict inequality)
        assert!(!conformal_causal_check(&p1, &p2));
    }
}
