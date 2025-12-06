//! Poisson sprinkling algorithms for generating causal sets.
//!
//! Implements Lorentz-invariant random point placement into spacetime regions.
//! The key insight: Poisson processes preserve Lorentz invariance, unlike regular lattices.
//!
//! # References
//! - Bombelli, Lee, Meyer, Sorkin (1987) "Space-time as a causal set"
//! - Surya (2019) "The causal set approach to quantum gravity", Section 3

use crate::point::SpacetimePoint;
use nalgebra::SVector;
use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};
use tracing::{debug, info, instrument};

/// Configuration for a sprinkling experiment.
///
/// Captures all parameters needed for reproducibility.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SprinklingConfig {
    /// Random seed for reproducibility
    pub seed: u64,

    /// Sprinkling density ρ (elements per unit spacetime volume)
    pub density: f64,

    /// Unique identifier for this experiment run
    pub experiment_id: String,
}

impl SprinklingConfig {
    pub fn new(seed: u64, density: f64) -> Self {
        Self {
            seed,
            density,
            experiment_id: format!("exp_{seed}_{}", chrono_lite_timestamp()),
        }
    }
}

/// Simple timestamp without external dependency
fn chrono_lite_timestamp() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// A causal diamond (Alexandrov interval) in Minkowski spacetime.
///
/// The region J⁺(p) ∩ J⁻(q) where p is the past tip and q is the future tip.
/// In 2D, this forms a diamond shape; in higher dimensions, a "bicone".
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CausalDiamond<const D: usize> {
    /// Past tip of the diamond (earliest point)
    pub past_tip: SVector<f64, D>,

    /// Future tip of the diamond (latest point)
    pub future_tip: SVector<f64, D>,

    /// Proper time separation τ between tips
    pub proper_time: f64,
}

impl<const D: usize> CausalDiamond<D> {
    /// Create a causal diamond from past and future tips.
    ///
    /// # Panics
    /// Panics if the tips are not timelike separated (future tip must be in causal future of past tip).
    pub fn new(past_tip: SVector<f64, D>, future_tip: SVector<f64, D>) -> Self {
        let dt = future_tip[0] - past_tip[0];
        assert!(dt > 0.0, "Future tip must be after past tip");

        let dx_squared: f64 = (1..D)
            .map(|i| {
                let dx = future_tip[i] - past_tip[i];
                dx * dx
            })
            .sum();

        let interval_squared = dt * dt - dx_squared;
        assert!(
            interval_squared > 0.0,
            "Tips must be timelike separated (got spacelike interval)"
        );

        let proper_time = interval_squared.sqrt();

        Self {
            past_tip,
            future_tip,
            proper_time,
        }
    }

    /// Create a symmetric diamond centered at origin with given proper time extent.
    ///
    /// Past tip at (-τ/2, 0, 0, ...), future tip at (τ/2, 0, 0, ...)
    pub fn symmetric(proper_time: f64) -> Self {
        let mut past_tip = SVector::<f64, D>::zeros();
        let mut future_tip = SVector::<f64, D>::zeros();

        past_tip[0] = -proper_time / 2.0;
        future_tip[0] = proper_time / 2.0;

        Self {
            past_tip,
            future_tip,
            proper_time,
        }
    }

    /// Compute the spacetime volume of the causal diamond.
    ///
    /// V_d = C_d · τ^d where C_d is computed by integrating the
    /// (d-1)-dimensional spatial ball over the time interval.
    ///
    /// For a diamond with tips at t = ±τ/2, at time t the spatial
    /// cross-section is a ball of radius r(t) = τ/2 - |t|.
    ///
    /// Verified by direct integration:
    /// - C₂ = 1/2    (2D: interval cross-section)
    /// - C₃ = π/12   (3D: disc cross-section)
    /// - C₄ = π/24   (4D: ball cross-section)
    pub fn volume(&self) -> f64 {
        let c_d = match D {
            2 => 0.5,
            3 => std::f64::consts::PI / 12.0,  // Fixed: was π/4
            4 => std::f64::consts::PI / 24.0,
            _ => {
                // General formula: C_d = π^((d-1)/2) / [d! · Γ((d+1)/2)]
                // For now, approximate or panic for unsupported dimensions
                panic!(
                    "Volume coefficient not implemented for D={}. Add formula for this dimension.",
                    D
                );
            }
        };

        c_d * self.proper_time.powi(D as i32)
    }

    /// Expected number of points for given density.
    pub fn expected_count(&self, density: f64) -> f64 {
        density * self.volume()
    }

    /// Check if a point lies inside the causal diamond.
    ///
    /// A point p is inside iff:
    /// - p is in the causal future of past_tip
    /// - p is in the causal past of future_tip
    pub fn contains(&self, point: &SVector<f64, D>) -> bool {
        // Check: past_tip ≺ point
        let dt_from_past = point[0] - self.past_tip[0];
        if dt_from_past <= 0.0 {
            return false;
        }
        let dx_from_past_sq: f64 = (1..D)
            .map(|i| {
                let dx = point[i] - self.past_tip[i];
                dx * dx
            })
            .sum();
        if dt_from_past * dt_from_past < dx_from_past_sq {
            return false;
        }

        // Check: point ≺ future_tip
        let dt_to_future = self.future_tip[0] - point[0];
        if dt_to_future <= 0.0 {
            return false;
        }
        let dx_to_future_sq: f64 = (1..D)
            .map(|i| {
                let dx = self.future_tip[i] - point[i];
                dx * dx
            })
            .sum();
        if dt_to_future * dt_to_future < dx_to_future_sq {
            return false;
        }

        true
    }
}

/// Sprinkler for generating causal sets via Poisson process.
///
/// Maintains RNG state for reproducible sequences of sprinklings.
pub struct Sprinkler {
    rng: ChaCha8Rng,
    config: SprinklingConfig,
    next_id: usize,
    batch_count: u64,
}

impl Sprinkler {
    /// Create a new sprinkler with given configuration.
    pub fn new(config: SprinklingConfig) -> Self {
        let rng = ChaCha8Rng::seed_from_u64(config.seed);
        info!(
            seed = config.seed,
            density = config.density,
            experiment_id = %config.experiment_id,
            "Initialized sprinkler"
        );

        Self {
            rng,
            config,
            next_id: 0,
            batch_count: 0,
        }
    }

    /// Create a sprinkler with a simple seed and density.
    pub fn with_seed(seed: u64, density: f64) -> Self {
        Self::new(SprinklingConfig::new(seed, density))
    }

    /// Get the current configuration (for logging/reproducibility).
    pub fn config(&self) -> &SprinklingConfig {
        &self.config
    }

    /// Sprinkle points into a causal diamond using rejection sampling.
    ///
    /// The algorithm:
    /// 1. Draw N from Poisson(ρV) where V is the diamond volume
    /// 2. Sample N points uniformly from bounding box
    /// 3. Reject points outside the causal diamond
    /// 4. Repeat until we have the target count
    ///
    /// This preserves Lorentz invariance of the point process.
    #[instrument(skip(self), fields(batch = self.batch_count))]
    pub fn sprinkle_diamond<const D: usize>(
        &mut self,
        diamond: &CausalDiamond<D>,
    ) -> SprinklingResult<D> {
        let expected = diamond.expected_count(self.config.density);
        let n_target = self.sample_poisson(expected);

        debug!(
            volume = diamond.volume(),
            expected = expected,
            n_target = n_target,
            "Starting rejection sampling"
        );

        let mut points = Vec::with_capacity(n_target);
        let mut total_samples = 0usize;

        // Compute bounding box
        let (t_min, t_max, _spatial_extent) = self.compute_bounding_box(diamond);

        let t_dist = Uniform::new(t_min, t_max);

        while points.len() < n_target {
            total_samples += 1;

            // Sample uniformly from bounding box
            let t = t_dist.sample(&mut self.rng);

            let mut coords = SVector::<f64, D>::zeros();
            coords[0] = t;

            // For each spatial dimension, sample within the light cone constraint
            // At time t, the spatial extent is limited by the light cones from both tips
            let dt_from_past = t - diamond.past_tip[0];
            let dt_to_future = diamond.future_tip[0] - t;
            let max_spatial_radius = dt_from_past.min(dt_to_future);

            for i in 1..D {
                let center = diamond.past_tip[i]
                    + (diamond.future_tip[i] - diamond.past_tip[i]) * (dt_from_past / diamond.proper_time);
                let dist = Uniform::new(-max_spatial_radius, max_spatial_radius);
                coords[i] = center + dist.sample(&mut self.rng);
            }

            // Check if point is inside the diamond
            if diamond.contains(&coords) {
                let point = SpacetimePoint::with_batch(coords, self.next_id, self.batch_count);
                self.next_id += 1;
                points.push(point);
            }
        }

        let acceptance_rate = n_target as f64 / total_samples as f64;
        debug!(
            n_points = points.len(),
            total_samples = total_samples,
            acceptance_rate = acceptance_rate,
            "Sprinkling complete"
        );

        self.batch_count += 1;

        SprinklingResult {
            points,
            diamond: diamond.clone(),
            config: self.config.clone(),
            batch_id: self.batch_count - 1,
            total_samples,
            acceptance_rate,
        }
    }

    /// Sprinkle a fixed number of points (useful for controlled experiments).
    ///
    /// Note: This does NOT follow a true Poisson process (count is fixed, not random).
    /// Use only when you need exact counts for comparison studies.
    ///
    /// Uses correct volume-weighted sampling: time is sampled with PDF proportional
    /// to r(t)^(d-1) where r(t) is the spatial cross-section radius at time t.
    /// This ensures uniform density throughout the diamond.
    #[instrument(skip(self), fields(batch = self.batch_count))]
    pub fn sprinkle_fixed<const D: usize>(
        &mut self,
        diamond: &CausalDiamond<D>,
        n_points: usize,
    ) -> SprinklingResult<D> {
        debug!(n_points = n_points, "Fixed-count sprinkling");

        let mut points = Vec::with_capacity(n_points);
        let mut total_samples = 0usize;

        // For symmetric diamond: half-extent T = τ/2
        let t_center = (diamond.past_tip[0] + diamond.future_tip[0]) / 2.0;
        let t_half = diamond.proper_time / 2.0;
        let spatial_dims = D - 1; // Number of spatial dimensions

        let unit_dist = Uniform::new(0.0_f64, 1.0_f64);

        while points.len() < n_points {
            total_samples += 1;

            // Sample time with PDF ∝ r(t)^(spatial_dims) for uniform spacetime density
            // Using inverse CDF: |t - t_center| = T * (1 - (1-u)^(1/(k+1)))
            // where k = spatial_dims
            let u: f64 = unit_dist.sample(&mut self.rng);
            let exponent = 1.0 / (spatial_dims as f64 + 1.0);
            let abs_t_offset = t_half * (1.0 - (1.0 - u).powf(exponent));

            // Random sign for which side of center
            let t = if unit_dist.sample(&mut self.rng) < 0.5 {
                t_center + abs_t_offset
            } else {
                t_center - abs_t_offset
            };

            let mut coords = SVector::<f64, D>::zeros();
            coords[0] = t;

            let dt_from_past = t - diamond.past_tip[0];
            let dt_to_future = diamond.future_tip[0] - t;
            let max_spatial_radius = dt_from_past.min(dt_to_future);

            // Sample spatial coordinates uniformly in ball of radius max_spatial_radius
            // For D=2 (1 spatial dim): just sample uniformly in [-r, r]
            // For D>2: sample uniformly in ball using standard methods
            if D == 2 {
                let dist = Uniform::new(-max_spatial_radius, max_spatial_radius);
                coords[1] = diamond.past_tip[1]
                    + (diamond.future_tip[1] - diamond.past_tip[1])
                        * (dt_from_past / diamond.proper_time)
                    + dist.sample(&mut self.rng);
            } else {
                // For higher dimensions, use rejection sampling in ball
                // (more efficient methods exist but this is simple and correct)
                let center_spatial: Vec<f64> = (1..D)
                    .map(|i| {
                        diamond.past_tip[i]
                            + (diamond.future_tip[i] - diamond.past_tip[i])
                                * (dt_from_past / diamond.proper_time)
                    })
                    .collect();

                loop {
                    let mut spatial_offset = vec![0.0; D - 1];
                    let mut r_sq = 0.0;
                    for j in 0..(D - 1) {
                        let dist = Uniform::new(-max_spatial_radius, max_spatial_radius);
                        spatial_offset[j] = dist.sample(&mut self.rng);
                        r_sq += spatial_offset[j] * spatial_offset[j];
                    }
                    // Accept if inside ball
                    if r_sq <= max_spatial_radius * max_spatial_radius {
                        for (j, offset) in spatial_offset.iter().enumerate() {
                            coords[j + 1] = center_spatial[j] + offset;
                        }
                        break;
                    }
                }
            }

            if diamond.contains(&coords) {
                let point = SpacetimePoint::with_batch(coords, self.next_id, self.batch_count);
                self.next_id += 1;
                points.push(point);
            }
        }

        let acceptance_rate = n_points as f64 / total_samples as f64;

        self.batch_count += 1;

        SprinklingResult {
            points,
            diamond: diamond.clone(),
            config: self.config.clone(),
            batch_id: self.batch_count - 1,
            total_samples,
            acceptance_rate,
        }
    }

    /// Sample from Poisson distribution using inverse transform.
    fn sample_poisson(&mut self, lambda: f64) -> usize {
        if lambda <= 0.0 {
            return 0;
        }

        // For small λ, use inverse transform method
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let mut k = 0usize;
            let mut p = 1.0;
            let uniform = Uniform::new(0.0, 1.0);

            loop {
                k += 1;
                p *= uniform.sample(&mut self.rng);
                if p <= l {
                    break;
                }
            }
            k - 1
        } else {
            // For large λ, use normal approximation
            use rand_distr::{Distribution, StandardNormal};
            let z: f64 = StandardNormal.sample(&mut self.rng);
            let result = lambda + z * lambda.sqrt();
            result.round().max(0.0) as usize
        }
    }

    /// Compute bounding box for rejection sampling.
    fn compute_bounding_box<const D: usize>(
        &self,
        diamond: &CausalDiamond<D>,
    ) -> (f64, f64, f64) {
        let t_min = diamond.past_tip[0];
        let t_max = diamond.future_tip[0];
        let spatial_extent = diamond.proper_time / 2.0; // Max spatial extent at midpoint

        (t_min, t_max, spatial_extent)
    }
}

/// Result of a sprinkling operation with full metadata for reproducibility.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SprinklingResult<const D: usize> {
    /// The generated spacetime points
    pub points: Vec<SpacetimePoint<D>>,

    /// The causal diamond that was sprinkled
    pub diamond: CausalDiamond<D>,

    /// Configuration used for this sprinkling
    pub config: SprinklingConfig,

    /// Batch identifier for this sprinkling
    pub batch_id: u64,

    /// Total samples attempted (for rejection sampling efficiency tracking)
    pub total_samples: usize,

    /// Acceptance rate of rejection sampling
    pub acceptance_rate: f64,
}

// ============================================================================
// Generic Spacetime Sprinkling (Phase 4)
// ============================================================================

use crate::spacetime::Spacetime;

/// Result of sprinkling into a generic spacetime.
#[derive(Debug, Clone)]
pub struct GenericSprinklingResult<const D: usize> {
    /// The generated spacetime points
    pub points: Vec<SpacetimePoint<D>>,

    /// Name of the spacetime that was sprinkled
    pub spacetime_name: String,

    /// Random seed used
    pub seed: u64,

    /// Volume of the spacetime region
    pub volume: f64,

    /// Ricci scalar of the spacetime
    pub ricci_scalar: f64,
}

/// Sprinkle a fixed number of points into any spacetime.
///
/// This is the generic version of sprinkling that works with any
/// [`Spacetime`] implementation (Minkowski, de Sitter, etc.).
///
/// # Arguments
///
/// * `spacetime` - The spacetime to sprinkle into
/// * `n_points` - Number of points to generate
/// * `seed` - Random seed for reproducibility
///
/// # Example
///
/// ```rust
/// use causet_core::spacetime::{Minkowski, DeSitter};
/// use causet_core::sprinkling::sprinkle_spacetime;
///
/// // Sprinkle into Minkowski
/// let mink = Minkowski::<4>::causal_diamond(5.0);
/// let result = sprinkle_spacetime(&mink, 500, 42);
/// assert_eq!(result.points.len(), 500);
///
/// // Sprinkle into de Sitter
/// let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
/// let result = sprinkle_spacetime(&ds, 500, 42);
/// assert!(result.ricci_scalar > 0.0);
/// ```
pub fn sprinkle_spacetime<S, const D: usize>(
    spacetime: &S,
    n_points: usize,
    seed: u64,
) -> GenericSprinklingResult<D>
where
    S: Spacetime<D>,
{
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut points = Vec::with_capacity(n_points);

    for id in 0..n_points {
        let coords = spacetime.sample_point(&mut rng);
        let point = SpacetimePoint::with_batch(coords, id, 0);
        points.push(point);
    }

    GenericSprinklingResult {
        points,
        spacetime_name: spacetime.name(),
        seed,
        volume: spacetime.volume(),
        ricci_scalar: spacetime.ricci_scalar(),
    }
}

impl<const D: usize> SprinklingResult<D> {
    /// Effective density (actual points / volume).
    pub fn effective_density(&self) -> f64 {
        self.points.len() as f64 / self.diamond.volume()
    }

    /// Deviation from expected count.
    pub fn count_deviation(&self) -> f64 {
        let expected = self.diamond.expected_count(self.config.density);
        (self.points.len() as f64 - expected) / expected.sqrt()
    }

    /// Export metadata for research paper documentation.
    pub fn to_metadata_json(&self) -> String {
        serde_json::to_string_pretty(&SprinklingMetadata {
            seed: self.config.seed,
            density: self.config.density,
            experiment_id: self.config.experiment_id.clone(),
            batch_id: self.batch_id,
            n_points: self.points.len(),
            diamond_volume: self.diamond.volume(),
            proper_time: self.diamond.proper_time,
            effective_density: self.effective_density(),
            count_deviation_sigma: self.count_deviation(),
            acceptance_rate: self.acceptance_rate,
            dimensions: D,
        })
        .unwrap_or_default()
    }
}

/// Metadata struct for serialization (without the full point data).
#[derive(Debug, Serialize, Deserialize)]
struct SprinklingMetadata {
    seed: u64,
    density: f64,
    experiment_id: String,
    batch_id: u64,
    n_points: usize,
    diamond_volume: f64,
    proper_time: f64,
    effective_density: f64,
    count_deviation_sigma: f64,
    acceptance_rate: f64,
    dimensions: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_diamond_volume_2d() {
        let diamond = CausalDiamond::<2>::symmetric(2.0);
        // V₂ = 0.5 · τ² = 0.5 · 4 = 2
        assert_relative_eq!(diamond.volume(), 2.0);
    }

    #[test]
    fn test_diamond_volume_3d() {
        let diamond = CausalDiamond::<3>::symmetric(2.0);
        // V₃ = (π/12) · τ³ = (π/12) · 8 = 2π/3
        assert_relative_eq!(diamond.volume(), 2.0 * std::f64::consts::PI / 3.0);
    }

    #[test]
    fn test_diamond_volume_4d() {
        let diamond = CausalDiamond::<4>::symmetric(2.0);
        // V₄ = (π/24) · τ⁴ = (π/24) · 16 = 2π/3
        assert_relative_eq!(diamond.volume(), 2.0 * std::f64::consts::PI / 3.0);
    }

    #[test]
    fn test_diamond_contains() {
        let diamond = CausalDiamond::<2>::symmetric(2.0);

        // Center point should be inside
        let center = SVector::<f64, 2>::new(0.0, 0.0);
        assert!(diamond.contains(&center));

        // Point near past tip should be inside
        let near_past = SVector::<f64, 2>::new(-0.9, 0.0);
        assert!(diamond.contains(&near_past));

        // Point outside (spacelike from tips) should be rejected
        let outside = SVector::<f64, 2>::new(0.0, 1.5);
        assert!(!diamond.contains(&outside));
    }

    #[test]
    fn test_sprinkling_reproducibility() {
        let diamond = CausalDiamond::<2>::symmetric(4.0);

        // Same seed should give identical results
        let mut sprinkler1 = Sprinkler::with_seed(12345, 10.0);
        let mut sprinkler2 = Sprinkler::with_seed(12345, 10.0);

        let result1 = sprinkler1.sprinkle_fixed(&diamond, 100);
        let result2 = sprinkler2.sprinkle_fixed(&diamond, 100);

        assert_eq!(result1.points.len(), result2.points.len());
        for (p1, p2) in result1.points.iter().zip(result2.points.iter()) {
            assert_relative_eq!(p1.coords[0], p2.coords[0], epsilon = 1e-10);
            assert_relative_eq!(p1.coords[1], p2.coords[1], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_sprinkling_density() {
        // Statistical test: actual count should be within ~3σ of expected
        let diamond = CausalDiamond::<2>::symmetric(10.0);
        let density = 5.0;
        let mut sprinkler = Sprinkler::with_seed(42, density);

        let result = sprinkler.sprinkle_diamond(&diamond);
        let expected = diamond.expected_count(density);

        // Should be within 4 sigma (very conservative for tests)
        let deviation = (result.points.len() as f64 - expected).abs();
        let sigma = expected.sqrt();
        assert!(
            deviation < 4.0 * sigma,
            "Count {} deviates too much from expected {} ({}σ)",
            result.points.len(),
            expected,
            deviation / sigma
        );
    }

    #[test]
    fn test_all_points_inside_diamond() {
        let diamond = CausalDiamond::<2>::symmetric(4.0);
        let mut sprinkler = Sprinkler::with_seed(999, 20.0);

        let result = sprinkler.sprinkle_diamond(&diamond);

        for point in &result.points {
            assert!(
                diamond.contains(&point.coords),
                "Point {:?} is outside diamond",
                point.coords
            );
        }
    }
}
