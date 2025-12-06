//! Flat Minkowski spacetime implementation.
//!
//! Minkowski spacetime is the flat spacetime of special relativity with metric:
//!
//! ```text
//! ds² = -dt² + dx² + dy² + dz² + ...
//! ```
//!
//! Sprinkling is performed in a causal diamond (Alexandrov interval) centered
//! at the origin with proper time extent τ.

use super::{conformal_causal_check, Spacetime};
use nalgebra::SVector;
use rand::distributions::{Distribution, Uniform};
use rand::Rng;

/// Flat Minkowski spacetime with sprinkling in a causal diamond.
///
/// The causal diamond is centered at the origin with time coordinate
/// ranging from -τ/2 to +τ/2, where τ is the proper time (half_height × 2).
///
/// # Volume Formula
///
/// The volume of a causal diamond in d dimensions is V_d = C_d × τ^d where:
/// - C₂ = 1/2
/// - C₃ = π/12
/// - C₄ = π/24
///
/// # Example
///
/// ```rust
/// use causet_core::spacetime::{Minkowski, Spacetime};
///
/// // Create a 4D Minkowski diamond with proper time τ = 10
/// let mink = Minkowski::<4>::causal_diamond(5.0); // half_height = τ/2
/// assert_eq!(mink.ricci_scalar(), 0.0); // Flat spacetime
/// ```
#[derive(Clone, Debug)]
pub struct Minkowski<const D: usize> {
    /// Half the proper time extent (τ/2)
    half_height: f64,
}

impl<const D: usize> Minkowski<D> {
    /// Create a Minkowski causal diamond with given half-height.
    ///
    /// The diamond spans time coordinate from -half_height to +half_height,
    /// giving a total proper time of 2 × half_height.
    ///
    /// # Arguments
    ///
    /// * `half_height` - Half the proper time extent (τ/2)
    ///
    /// # Panics
    ///
    /// Panics if half_height is not positive.
    pub fn causal_diamond(half_height: f64) -> Self {
        assert!(half_height > 0.0, "half_height must be positive");
        Self { half_height }
    }

    /// Get the half-height (τ/2) of the causal diamond.
    pub fn half_height(&self) -> f64 {
        self.half_height
    }

    /// Get the full proper time extent τ.
    pub fn proper_time(&self) -> f64 {
        2.0 * self.half_height
    }

    /// Volume coefficient C_d for dimension D.
    fn volume_coefficient() -> f64 {
        match D {
            2 => 0.5,
            3 => std::f64::consts::PI / 12.0,
            4 => std::f64::consts::PI / 24.0,
            _ => panic!(
                "Volume coefficient not implemented for D={}. \
                 General formula: C_d = π^((d-1)/2) / [d! × Γ((d+1)/2)]",
                D
            ),
        }
    }
}

impl<const D: usize> Spacetime<D> for Minkowski<D> {
    fn name(&self) -> String {
        format!("Minkowski-{}D", D)
    }

    fn sample_point(&self, rng: &mut impl Rng) -> SVector<f64, D> {
        let t_half = self.half_height;
        let spatial_dims = D - 1;
        let unit_dist = Uniform::new(0.0_f64, 1.0_f64);

        // Sample time with PDF ∝ r(t)^(spatial_dims) for uniform spacetime density
        // At time t, the spatial cross-section radius is r(t) = half_height - |t|
        // Using inverse CDF: |t| = T × (1 - (1-u)^(1/(k+1))) where k = spatial_dims
        let u: f64 = unit_dist.sample(rng);
        let exponent = 1.0 / (spatial_dims as f64 + 1.0);
        let abs_t = t_half * (1.0 - (1.0 - u).powf(exponent));

        // Random sign for which half of the diamond
        let t = if rng.gen::<bool>() {
            abs_t
        } else {
            -abs_t
        };

        let mut coords = SVector::<f64, D>::zeros();
        coords[0] = t;

        // Maximum spatial radius at this time
        let max_r = t_half - abs_t;

        // Sample spatial coordinates uniformly in ball of radius max_r
        if D == 2 {
            // 1D spatial: uniform in [-r, r]
            coords[1] = rng.gen_range(-max_r..max_r);
        } else {
            // Higher dimensions: rejection sampling in ball
            loop {
                let mut r_sq = 0.0;
                for i in 1..D {
                    let x = rng.gen_range(-max_r..max_r);
                    coords[i] = x;
                    r_sq += x * x;
                }
                if r_sq <= max_r * max_r {
                    break;
                }
            }
        }

        coords
    }

    fn causally_precedes(&self, p1: &SVector<f64, D>, p2: &SVector<f64, D>) -> bool {
        conformal_causal_check(p1, p2)
    }

    fn ricci_scalar(&self) -> f64 {
        0.0 // Flat spacetime
    }

    fn volume(&self) -> f64 {
        let tau = self.proper_time();
        Self::volume_coefficient() * tau.powi(D as i32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_minkowski_volume_2d() {
        let mink = Minkowski::<2>::causal_diamond(1.0); // τ = 2
        // V = 0.5 × 2² = 2
        assert_relative_eq!(mink.volume(), 2.0);
    }

    #[test]
    fn test_minkowski_volume_4d() {
        let mink = Minkowski::<4>::causal_diamond(1.0); // τ = 2
        // V = (π/24) × 2⁴ = (π/24) × 16 = 2π/3
        assert_relative_eq!(mink.volume(), 2.0 * std::f64::consts::PI / 3.0);
    }

    #[test]
    fn test_minkowski_ricci_scalar() {
        let mink = Minkowski::<4>::causal_diamond(5.0);
        assert_eq!(mink.ricci_scalar(), 0.0);
    }

    #[test]
    fn test_minkowski_name() {
        let mink = Minkowski::<4>::causal_diamond(5.0);
        assert_eq!(mink.name(), "Minkowski-4D");
    }

    #[test]
    fn test_sampled_points_are_timelike_from_origin() {
        let mink = Minkowski::<4>::causal_diamond(5.0);
        let mut rng = ChaCha8Rng::seed_from_u64(42);

        // Sample many points and verify they're inside the light cone from tips
        for _ in 0..100 {
            let p = mink.sample_point(&mut rng);
            let t = p[0];
            let r_sq: f64 = (1..4).map(|i| p[i] * p[i]).sum();
            let r = r_sq.sqrt();

            // Point should be inside diamond: |t| + r < half_height
            assert!(
                t.abs() + r <= mink.half_height() + 1e-10,
                "Point {:?} outside diamond",
                p
            );
        }
    }

    #[test]
    fn test_causal_relation() {
        let mink = Minkowski::<4>::causal_diamond(5.0);

        // Timelike separated (causal)
        let p1 = SVector::<f64, 4>::new(0.0, 0.0, 0.0, 0.0);
        let p2 = SVector::<f64, 4>::new(2.0, 1.0, 0.0, 0.0);
        assert!(mink.causally_precedes(&p1, &p2));

        // Spacelike separated (not causal)
        let p3 = SVector::<f64, 4>::new(1.0, 2.0, 0.0, 0.0);
        assert!(!mink.causally_precedes(&p1, &p3));
    }
}
