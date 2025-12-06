//! de Sitter spacetime implementation in conformal coordinates.
//!
//! de Sitter spacetime is a maximally symmetric spacetime with positive
//! constant curvature. It describes an exponentially expanding universe
//! and is relevant for inflationary cosmology.
//!
//! # Conformal Coordinates
//!
//! We use conformal (flat slicing) coordinates where the metric is:
//!
//! ```text
//! ds² = (1/H²η²)(-dη² + dx⃗²)
//! ```
//!
//! where:
//! - η is conformal time, ranging from -∞ (past) to 0⁻ (future infinity)
//! - H is the Hubble parameter (sets the curvature scale)
//! - x⃗ are comoving spatial coordinates
//!
//! # Key Properties
//!
//! - **Conformally flat**: Causal structure is identical to Minkowski
//! - **Volume element**: √(-g) = 1/(H|η|)^d
//! - **Ricci scalar**: R = d(d-1)H² (constant, positive)

use super::{conformal_causal_check, Spacetime};
use nalgebra::SVector;
use rand::Rng;

/// de Sitter spacetime in conformal coordinates.
///
/// The conformal metric is ds² = (1/H²η²)(-dη² + dx⃗²) where η < 0.
///
/// # Physics
///
/// de Sitter spacetime has:
/// - Constant positive scalar curvature R = d(d-1)H²
/// - Exponential expansion (in proper time)
/// - Conformal structure identical to Minkowski
///
/// # Sprinkling Region
///
/// Points are sampled in the region:
/// - Conformal time: η ∈ [eta_min, eta_max] (both negative)
/// - Spatial: x_i ∈ [-L, +L] for each spatial dimension
///
/// # Example
///
/// ```rust
/// use causet_core::spacetime::{DeSitter, Spacetime};
///
/// // Create de Sitter with H = 0.1, conformal time from -10 to -1
/// let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
/// assert!(ds.ricci_scalar() > 0.0); // Positive curvature
/// ```
#[derive(Clone, Debug)]
pub struct DeSitter<const D: usize> {
    /// Hubble parameter (curvature scale)
    h: f64,
    /// Conformal time range (both negative: eta_min < eta_max < 0)
    eta_range: (f64, f64),
    /// Spatial half-width (spatial extent is [-L, +L] in each direction)
    spatial_half_width: f64,
}

impl<const D: usize> DeSitter<D> {
    /// Create a de Sitter spacetime patch with given parameters.
    ///
    /// # Arguments
    ///
    /// * `h` - Hubble parameter (positive, sets curvature: R = d(d-1)H²)
    /// * `eta_min` - Past conformal time bound (negative, e.g., -10.0)
    /// * `eta_max` - Future conformal time bound (negative, closer to 0, e.g., -1.0)
    /// * `spatial_half_width` - Spatial extent [-L, +L] in each direction
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - h ≤ 0
    /// - eta_min ≥ eta_max
    /// - eta_max ≥ 0
    /// - spatial_half_width ≤ 0
    pub fn new(h: f64, eta_min: f64, eta_max: f64, spatial_half_width: f64) -> Self {
        assert!(h > 0.0, "Hubble parameter must be positive, got {}", h);
        assert!(
            eta_min < eta_max,
            "eta_min must be less than eta_max: {} < {}",
            eta_min,
            eta_max
        );
        assert!(
            eta_max < 0.0,
            "Both eta bounds must be negative (conformal time), got eta_max = {}",
            eta_max
        );
        assert!(
            spatial_half_width > 0.0,
            "Spatial half-width must be positive"
        );

        Self {
            h,
            eta_range: (eta_min, eta_max),
            spatial_half_width,
        }
    }

    /// Get the Hubble parameter H.
    pub fn hubble_parameter(&self) -> f64 {
        self.h
    }

    /// Get the conformal time range (eta_min, eta_max).
    pub fn eta_range(&self) -> (f64, f64) {
        self.eta_range
    }

    /// Get the spatial half-width L.
    pub fn spatial_half_width(&self) -> f64 {
        self.spatial_half_width
    }

    /// Sample conformal time η with importance sampling.
    ///
    /// The volume element √(-g) = 1/(H|η|)^D, so we need to sample
    /// η with PDF ∝ |η|^(-D).
    ///
    /// Using inverse CDF method:
    /// - CDF: P(η) ∝ ∫ dη/|η|^D = |η|^(1-D)/(1-D)
    /// - Inverse: η = -(a + u(b-a))^(-1/(D-1))
    ///   where a = |η_min|^(1-D), b = |η_max|^(1-D)
    fn sample_eta(&self, u: f64) -> f64 {
        let (eta_min, eta_max) = self.eta_range;
        let exp = (D - 1) as i32;

        // For D=1, this doesn't make sense, but D≥2 for spacetime
        assert!(D >= 2, "Dimension must be at least 2");

        let a = eta_min.abs().powi(-exp);
        let b = eta_max.abs().powi(-exp);

        -(a + u * (b - a)).powf(-1.0 / exp as f64)
    }
}

impl<const D: usize> Spacetime<D> for DeSitter<D> {
    fn name(&self) -> String {
        format!("deSitter-{}D (H={:.3})", D, self.h)
    }

    fn sample_point(&self, rng: &mut impl Rng) -> SVector<f64, D> {
        // Sample conformal time with importance sampling
        let u: f64 = rng.gen();
        let eta = self.sample_eta(u);

        let mut coords = SVector::<f64, D>::zeros();
        coords[0] = eta;

        // Spatial coordinates: uniform in box (conformal coords are spatially flat)
        let l = self.spatial_half_width;
        for i in 1..D {
            coords[i] = rng.gen_range(-l..l);
        }

        coords
    }

    fn causally_precedes(&self, p1: &SVector<f64, D>, p2: &SVector<f64, D>) -> bool {
        // Conformal coordinates: causal structure identical to Minkowski
        conformal_causal_check(p1, p2)
    }

    fn ricci_scalar(&self) -> f64 {
        // de Sitter: R = d(d-1)H²
        // 4D: R = 12H²
        let d = D as f64;
        d * (d - 1.0) * self.h * self.h
    }

    fn volume(&self) -> f64 {
        // Proper volume: ∫ √(-g) d^D x = ∫ (1/(H|η|))^D dη d^(D-1)x
        // = (spatial volume) × ∫_{η_min}^{η_max} dη / (H|η|)^D
        let (eta_min, eta_max) = self.eta_range;
        let spatial_dims = D - 1;

        // Spatial volume: (2L)^(D-1)
        let spatial_vol = (2.0 * self.spatial_half_width).powi(spatial_dims as i32);

        // Time integral: ∫ dη/|η|^D
        // For η < 0: the antiderivative of |η|^(-D) is |η|^(1-D)/(D-1)
        // Evaluated: [|η_max|^(1-D) - |η_min|^(1-D)] / (D-1)
        // Since D > 1 and |η_max| < |η_min|, the numerator is positive
        let exp = 1 - D as i32;
        let time_integral =
            (eta_max.abs().powi(exp) - eta_min.abs().powi(exp)) / (D as f64 - 1.0);

        spatial_vol * time_integral / self.h.powi(D as i32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_de_sitter_ricci_scalar_4d() {
        let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
        // R = d(d-1)H² = 4×3×0.01 = 0.12 = 12H²
        assert_relative_eq!(ds.ricci_scalar(), 0.12);
    }

    #[test]
    fn test_de_sitter_ricci_scalar_2d() {
        let ds = DeSitter::<2>::new(0.5, -10.0, -1.0, 5.0);
        // R = d(d-1)H² = 2×1×0.25 = 0.5 = 2H²
        assert_relative_eq!(ds.ricci_scalar(), 0.5);
    }

    #[test]
    fn test_de_sitter_name() {
        let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
        assert!(ds.name().contains("deSitter"));
        assert!(ds.name().contains("4D"));
    }

    #[test]
    fn test_de_sitter_volume_positive() {
        let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
        assert!(ds.volume() > 0.0, "Volume should be positive");
    }

    #[test]
    fn test_de_sitter_volume_scales_with_h() {
        // Volume ∝ 1/H^D, so larger H → smaller volume
        let ds1 = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
        let ds2 = DeSitter::<4>::new(0.2, -10.0, -1.0, 5.0);

        // V ∝ 1/H^4, so V2/V1 = (H1/H2)^4 = (0.1/0.2)^4 = 1/16
        let ratio = ds2.volume() / ds1.volume();
        assert_relative_eq!(ratio, (0.1 / 0.2_f64).powi(4), epsilon = 1e-10);
    }

    #[test]
    fn test_sampled_eta_in_range() {
        let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
        let mut rng = ChaCha8Rng::seed_from_u64(42);

        for _ in 0..100 {
            let p = ds.sample_point(&mut rng);
            let eta = p[0];

            assert!(
                eta >= -10.0 && eta <= -1.0,
                "eta = {} outside range [-10, -1]",
                eta
            );
            assert!(eta < 0.0, "eta must be negative");
        }
    }

    #[test]
    fn test_sampled_spatial_in_range() {
        let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
        let mut rng = ChaCha8Rng::seed_from_u64(42);

        for _ in 0..100 {
            let p = ds.sample_point(&mut rng);

            for i in 1..4 {
                assert!(
                    p[i] >= -5.0 && p[i] <= 5.0,
                    "spatial coord {} = {} outside range [-5, 5]",
                    i,
                    p[i]
                );
            }
        }
    }

    #[test]
    fn test_causal_relation_same_as_minkowski() {
        let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);

        // Timelike separated (note: η increases toward 0, so -5 is "past" of -3)
        let p1 = SVector::<f64, 4>::new(-5.0, 0.0, 0.0, 0.0);
        let p2 = SVector::<f64, 4>::new(-3.0, 1.0, 0.0, 0.0); // Δη=2, |Δx|=1 → timelike
        assert!(ds.causally_precedes(&p1, &p2));

        // Spacelike separated
        let p3 = SVector::<f64, 4>::new(-4.0, 3.0, 0.0, 0.0); // Δη=1, |Δx|=3 → spacelike
        assert!(!ds.causally_precedes(&p1, &p3));
    }

    #[test]
    #[should_panic]
    fn test_de_sitter_invalid_h() {
        DeSitter::<4>::new(-0.1, -10.0, -1.0, 5.0);
    }

    #[test]
    #[should_panic]
    fn test_de_sitter_invalid_eta_order() {
        DeSitter::<4>::new(0.1, -1.0, -10.0, 5.0); // eta_min > eta_max
    }

    #[test]
    #[should_panic]
    fn test_de_sitter_positive_eta() {
        DeSitter::<4>::new(0.1, -10.0, 1.0, 5.0); // eta_max > 0
    }

    #[test]
    fn test_importance_sampling_distribution() {
        // Verify that importance sampling produces more points near η → 0
        // (where the volume element diverges)
        let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
        let mut rng = ChaCha8Rng::seed_from_u64(42);

        let n_samples = 1000;
        let mut near_zero_count = 0;
        let midpoint = -5.5;

        for _ in 0..n_samples {
            let p = ds.sample_point(&mut rng);
            if p[0] > midpoint {
                near_zero_count += 1;
            }
        }

        // Due to 1/|η|^4 weighting, most points should be near η = -1 (closer to 0)
        // The fraction near zero should be much larger than 0.5
        let fraction_near_zero = near_zero_count as f64 / n_samples as f64;
        assert!(
            fraction_near_zero > 0.8,
            "Expected most points near η=0, got fraction = {}",
            fraction_near_zero
        );
    }
}
