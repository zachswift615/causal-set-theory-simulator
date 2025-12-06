//! Phase 4 Validation Tests: Curved Spacetime Sprinkling
//!
//! These tests verify that the Benincasa-Dowker action detects curvature
//! in de Sitter spacetime while still giving ~0 for flat Minkowski.
//!
//! NOTE: The BD action has large fluctuations (σ ≈ 30×√N for small N).
//! Tests use many runs to get reliable statistics.

use causet_core::prelude::*;

/// Test 1: de Sitter gives positive BD action (most critical)
///
/// de Sitter has positive scalar curvature R = 12H², so the BD action
/// should be positive on average.
#[test]
fn de_sitter_positive_action() {
    let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);

    let n_runs = 20;
    let n_points = 300;
    let mut actions: Vec<f64> = Vec::with_capacity(n_runs);

    for seed in 0..n_runs {
        let result = sprinkle_spacetime(&ds, n_points, seed as u64);
        let causet = CausalSet::from_generic_sprinkling(result, &ds);
        actions.push(causet.bd_action_4d().action);
    }

    let mean: f64 = actions.iter().sum::<f64>() / n_runs as f64;
    let variance: f64 =
        actions.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / (n_runs - 1) as f64;
    let std_error = (variance / n_runs as f64).sqrt();

    println!(
        "de Sitter (H=0.1, N={}): mean S = {:.1}, SE = {:.1}",
        n_points, mean, std_error
    );

    // With large fluctuations, just check mean is positive
    // (The physics is correct if S > 0 for positive curvature)
    assert!(
        mean > 0.0,
        "de Sitter mean action should be positive, got mean={:.1}",
        mean
    );
}

/// Test 2: Different H values still give positive action
///
/// Note: In conformal coordinates, the causal structure is identical
/// regardless of H - only the volume element changes. With fixed N points,
/// the BD action doesn't directly scale with H (that would require fixed density).
/// Instead, we verify that de Sitter gives positive mean S for different H values.
#[test]
fn de_sitter_positive_for_various_h() {
    let h_values = [0.05, 0.15, 0.25];
    let n_points = 300;
    let n_runs = 15;

    for &h in &h_values {
        let ds = DeSitter::<4>::new(h, -10.0, -1.0, 5.0);
        let mut sum = 0.0;

        for seed in 0..n_runs {
            let result = sprinkle_spacetime(&ds, n_points, seed as u64);
            let causet = CausalSet::from_generic_sprinkling(result, &ds);
            sum += causet.bd_action_4d().action;
        }

        let mean = sum / n_runs as f64;
        let ricci = 12.0 * h * h;
        println!("H = {:.2}: mean S = {:.1}, R = {:.4}", h, mean, ricci);

        // Each H value should give positive mean S (positive curvature)
        assert!(
            mean > 0.0,
            "de Sitter (H={:.2}) should have positive mean action, got {:.1}",
            h, mean
        );
    }
}

/// Test 3: de Sitter mean action is positive (the key result)
///
/// This is the core validation: curved spacetime (de Sitter with R > 0)
/// should give a positive mean BD action, while Minkowski should give ~0.
#[test]
fn de_sitter_mean_positive_vs_minkowski() {
    let n_points = 300;
    let n_runs = 20;

    // Flat Minkowski
    let mink = Minkowski::<4>::causal_diamond(5.0);
    let mut mink_actions = Vec::new();
    for seed in 0..n_runs {
        let result = sprinkle_spacetime(&mink, n_points, seed as u64);
        let causet = CausalSet::from_generic_sprinkling(result, &mink);
        mink_actions.push(causet.bd_action_4d().action);
    }

    // Curved de Sitter
    let ds = DeSitter::<4>::new(0.1, -10.0, -1.0, 5.0);
    let mut ds_actions = Vec::new();
    for seed in 0..n_runs {
        let result = sprinkle_spacetime(&ds, n_points, seed as u64);
        let causet = CausalSet::from_generic_sprinkling(result, &ds);
        ds_actions.push(causet.bd_action_4d().action);
    }

    let mink_mean: f64 = mink_actions.iter().sum::<f64>() / n_runs as f64;
    let ds_mean: f64 = ds_actions.iter().sum::<f64>() / n_runs as f64;

    println!("Minkowski: mean S = {:.1}", mink_mean);
    println!("de Sitter: mean S = {:.1}", ds_mean);

    // Key physics: de Sitter (positive curvature) should have positive mean S
    // Minkowski mean should be ~0 (within fluctuations)
    assert!(
        ds_mean > 0.0,
        "de Sitter mean action should be positive, got {:.1}",
        ds_mean
    );

    // Also check that Minkowski is roughly centered around 0
    // (large fluctuations expected, so very loose bound)
    assert!(
        mink_mean.abs() < 500.0,
        "Minkowski mean action should be near zero, got {:.1}",
        mink_mean
    );
}

/// Test 4: Generic API gives same statistics as old API
///
/// The new generic sprinkling API with Minkowski should give
/// statistically equivalent results to the old CausalDiamond API.
#[test]
fn generic_api_matches_old_api_statistics() {
    let n_runs = 15;
    let n_points = 300;

    let mut old_fractions = Vec::new();
    let mut new_fractions = Vec::new();

    for seed in 0..n_runs {
        // Old API
        let diamond = CausalDiamond::<4>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(seed as u64, 1.0);
        let old_result = sprinkler.sprinkle_fixed(&diamond, n_points);
        let old_causet = CausalSet::from_sprinkling(old_result);
        old_fractions.push(old_causet.ordering_fraction());

        // New API
        let mink = Minkowski::<4>::causal_diamond(5.0);
        let new_result = sprinkle_spacetime(&mink, n_points, seed as u64);
        let new_causet = CausalSet::from_generic_sprinkling(new_result, &mink);
        new_fractions.push(new_causet.ordering_fraction());
    }

    let old_mean: f64 = old_fractions.iter().sum::<f64>() / n_runs as f64;
    let new_mean: f64 = new_fractions.iter().sum::<f64>() / n_runs as f64;

    println!("Old API mean ordering fraction: {:.4}", old_mean);
    println!("New API mean ordering fraction: {:.4}", new_mean);

    // Both should be close to 0.100 (theoretical for 4D)
    assert!(
        (old_mean - 0.100).abs() < 0.02,
        "Old API ordering fraction should be ~0.100, got {:.4}",
        old_mean
    );
    assert!(
        (new_mean - 0.100).abs() < 0.02,
        "New API ordering fraction should be ~0.100, got {:.4}",
        new_mean
    );
}

/// Test 5: 2D de Sitter produces valid causal sets
///
/// Sanity check that 2D de Sitter sprinkling produces non-trivial causal sets.
/// Note: The ordering fraction depends on region shape. Our rectangular region
/// in conformal coordinates gives f < 0.5 (the 0.5 value is for causal diamonds).
#[test]
fn de_sitter_2d_valid_causal_structure() {
    let ds = DeSitter::<2>::new(0.2, -10.0, -1.0, 5.0);

    let n_runs = 10;
    let n_points = 300;
    let mut fractions = Vec::new();

    for seed in 0..n_runs {
        let result = sprinkle_spacetime(&ds, n_points, seed as u64);
        let causet = CausalSet::from_generic_sprinkling(result, &ds);
        fractions.push(causet.ordering_fraction());
    }

    let mean: f64 = fractions.iter().sum::<f64>() / n_runs as f64;
    println!(
        "2D de Sitter (H=0.2, N={}): mean ordering fraction = {:.3}",
        n_points, mean
    );

    // Should be between 0 and 1 (valid causal structure)
    assert!(
        mean > 0.1 && mean < 0.9,
        "2D de Sitter ordering fraction should be non-trivial, got {:.3}",
        mean
    );
}
