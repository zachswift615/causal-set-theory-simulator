//! Phase 3: Benincasa-Dowker Action Validation
//!
//! This example demonstrates that the BD action computes correctly for
//! flat Minkowski spacetime, where we expect S ≈ 0 with O(√N) fluctuations.
//!
//! # Usage
//! ```bash
//! cargo run --example bd_action -p causet-core
//! ```

use causet_core::prelude::*;
use causet_core::sprinkling::{CausalDiamond, Sprinkler};

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║   Phase 3: Benincasa-Dowker Action for Flat Spacetime        ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("The BD action encodes scalar curvature. For flat Minkowski spacetime,");
    println!("we expect S ≈ 0 with statistical fluctuations of magnitude O(√N).");
    println!();

    println!("═══════════════════════════════════════════════════════════════");
    println!("2D Minkowski Spacetime: S^(2) = N - 2N₁ + 4N₂ - 2N₃");
    println!("═══════════════════════════════════════════════════════════════");
    println!();

    run_bd_validation::<2>("2D", 5, 500);

    println!();
    println!("═══════════════════════════════════════════════════════════════");
    println!("4D Minkowski Spacetime: S^(4) = N - N₁ + 9N₂ - 16N₃ + 8N₄");
    println!("═══════════════════════════════════════════════════════════════");
    println!();

    run_bd_validation::<4>("4D", 5, 500);

    println!();
    println!("═══════════════════════════════════════════════════════════════");
    println!("INTERPRETATION");
    println!("═══════════════════════════════════════════════════════════════");
    println!();
    println!("For flat spacetime:");
    println!("  • Raw action S fluctuates around 0");
    println!("  • Normalized action S/√N should be O(1)");
    println!("  • Fluctuations decrease as N increases");
    println!();
    println!("For curved spacetime (future work):");
    println!("  • S > 0 indicates positive curvature (de Sitter)");
    println!("  • S < 0 indicates negative curvature");
    println!("  • In continuum limit: S → ∫√(-g)R d⁴x (Einstein-Hilbert)");
}

fn run_bd_validation<const D: usize>(name: &str, n_trials: usize, n_points: usize) {
    let diamond = CausalDiamond::<D>::symmetric(10.0);

    println!(
        "Running {} trials with N={} points each\n",
        n_trials, n_points
    );
    println!("Trial    N       N₁      N₂      N₃      N₄      S        S/√N");
    println!("─────  ─────   ─────   ─────   ─────   ─────   ─────────  ──────");

    let mut actions = Vec::new();
    let mut normalized = Vec::new();

    for trial in 0..n_trials {
        let seed = 42 + trial as u64 * 1000;
        let mut sprinkler = Sprinkler::with_seed(seed, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, n_points);
        let causet = CausalSet::from_sprinkling(result);

        let action_result = if D == 2 {
            causet.bd_action_2d()
        } else {
            causet.bd_action_4d()
        };

        let counts = &action_result.interval_counts;
        let n1 = counts.n_k(0);
        let n2 = counts.n_k(1);
        let n3 = counts.n_k(2);
        let n4 = counts.n_k(3);

        println!(
            "{:>5}  {:>5}   {:>5}   {:>5}   {:>5}   {:>5}   {:>9.1}  {:>6.2}",
            trial + 1,
            causet.len(),
            n1,
            n2,
            n3,
            n4,
            action_result.action,
            action_result.normalized_action()
        );

        actions.push(action_result.action);
        normalized.push(action_result.normalized_action());
    }

    let mean_s: f64 = actions.iter().sum::<f64>() / actions.len() as f64;
    let mean_norm: f64 = normalized.iter().sum::<f64>() / normalized.len() as f64;
    let std_s: f64 = (actions
        .iter()
        .map(|&a| (a - mean_s).powi(2))
        .sum::<f64>()
        / (actions.len() - 1) as f64)
        .sqrt();
    let std_norm: f64 = (normalized
        .iter()
        .map(|&a| (a - mean_norm).powi(2))
        .sum::<f64>()
        / (normalized.len() - 1) as f64)
        .sqrt();

    println!();
    println!("Statistics for {} flat spacetime:", name);
    println!("  Mean S:     {:.1} ± {:.1}", mean_s, std_s);
    println!("  Mean S/√N:  {:.2} ± {:.2}", mean_norm, std_norm);
    println!(
        "  Expected:   S ≈ 0 with |S| < O(√N) ≈ {:.0}",
        (n_points as f64).sqrt() * 3.0
    );

    // Check if consistent with flat
    let expected_fluctuation = (n_points as f64).sqrt();
    if mean_s.abs() < 3.0 * expected_fluctuation {
        println!("  Result:     ✓ Consistent with flat spacetime (R = 0)");
    } else {
        println!(
            "  Result:     ⚠ Action larger than expected for flat spacetime"
        );
    }
}
