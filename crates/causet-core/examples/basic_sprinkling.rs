//! Basic Sprinkling Example — Phase 1 Validation
//!
//! This example demonstrates the core sprinkling engine by:
//! 1. Sprinkling points into 2D, 3D, and 4D causal diamonds
//! 2. Computing the ordering fraction for each dimension
//! 3. Reporting empirical observations for research documentation
//!
//! # Key Observable
//! The ordering fraction f = 2R/[N(N-1)] encodes spacetime dimension.
//! Higher dimensions yield lower ordering fractions.
//!
//! # Usage
//! ```bash
//! cargo run --example basic_sprinkling -p causet-core
//! cargo run --example basic_sprinkling -p causet-core -- --seed 12345
//! cargo run --example basic_sprinkling -p causet-core -- --n 500
//! ```

use causet_core::prelude::*;
use std::env;

fn main() {
    // Parse arguments
    let args: Vec<String> = env::args().collect();
    let seed = parse_arg(&args, "--seed").unwrap_or(42);
    let n_points = parse_arg(&args, "--n").unwrap_or(500) as usize;
    let n_trials = parse_arg(&args, "--trials").unwrap_or(5) as usize;

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║     CAUSET Phase 1: Basic Sprinkling Experiments            ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("Measuring ordering fractions in Minkowski spacetime diamonds.");
    println!();
    println!("Parameters:");
    println!("  • Seed: {} (--seed N)", seed);
    println!("  • Points per trial: {} (--n N)", n_points);
    println!("  • Trials per dimension: {} (--trials N)", n_trials);
    println!();

    // Run experiments for each dimension
    println!("═══════════════════════════════════════════════════════════════");
    let results_2d = run_experiment::<2>("2D (1+1)", seed, n_points, n_trials);

    println!();
    println!("═══════════════════════════════════════════════════════════════");
    let results_3d = run_experiment::<3>("3D (2+1)", seed + 1000, n_points, n_trials);

    println!();
    println!("═══════════════════════════════════════════════════════════════");
    let results_4d = run_experiment::<4>("4D (3+1)", seed + 2000, n_points, n_trials);

    // Summary comparison
    println!();
    println!("═══════════════════════════════════════════════════════════════");
    println!("SUMMARY: Dimension vs Ordering Fraction");
    println!("═══════════════════════════════════════════════════════════════");
    println!();
    println!("Dim   f_observed ± σ      Links/N");
    println!("───   ──────────────      ───────");
    print_summary("2D", &results_2d);
    print_summary("3D", &results_3d);
    print_summary("4D", &results_4d);

    println!();
    println!("Key observation: Ordering fraction DECREASES with dimension.");
    println!("This is the Myrheim-Meyer effect: higher dimensions → sparser causal structure.");
    println!();

    // Output JSON for research documentation
    println!("═══════════════════════════════════════════════════════════════");
    println!("Research Data (JSON for reproducibility):");
    println!("═══════════════════════════════════════════════════════════════");

    println!("{{");
    println!("  \"experiment\": \"basic_sprinkling_phase1\",");
    println!("  \"parameters\": {{");
    println!("    \"seed\": {},", seed);
    println!("    \"n_points\": {},", n_points);
    println!("    \"n_trials\": {}", n_trials);
    println!("  }},");
    println!("  \"results\": {{");
    println!("    \"2D\": {{ \"f_mean\": {:.6}, \"f_std\": {:.6}, \"links_per_element\": {:.3} }},",
             results_2d.0, results_2d.1, results_2d.2);
    println!("    \"3D\": {{ \"f_mean\": {:.6}, \"f_std\": {:.6}, \"links_per_element\": {:.3} }},",
             results_3d.0, results_3d.1, results_3d.2);
    println!("    \"4D\": {{ \"f_mean\": {:.6}, \"f_std\": {:.6}, \"links_per_element\": {:.3} }}",
             results_4d.0, results_4d.1, results_4d.2);
    println!("  }}");
    println!("}}");
}

fn parse_arg(args: &[String], flag: &str) -> Option<u64> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1))
        .and_then(|s| s.parse().ok())
}

/// Returns (mean_f, std_f, mean_link_ratio)
fn run_experiment<const D: usize>(
    name: &str,
    base_seed: u64,
    n_points: usize,
    n_trials: usize,
) -> (f64, f64, f64) {
    println!("{} Minkowski Spacetime", name);
    println!();

    let proper_time = 10.0;
    let diamond = CausalDiamond::<D>::symmetric(proper_time);

    println!("Diamond proper time: {:.1}", proper_time);
    println!("Diamond volume: {:.4}", diamond.volume());
    println!();

    let mut ordering_fractions = Vec::with_capacity(n_trials);
    let mut link_ratios = Vec::with_capacity(n_trials);

    println!("Trial  Points  Relations  Links    f_obs    Links/N");
    println!("─────  ──────  ─────────  ─────    ─────    ───────");

    for trial in 0..n_trials {
        let seed = base_seed + trial as u64;
        let mut sprinkler = Sprinkler::with_seed(seed, 1.0);
        let result = sprinkler.sprinkle_fixed(&diamond, n_points);
        let causet = CausalSet::from_sprinkling(result);

        let f_obs = causet.ordering_fraction();
        let n = causet.len();
        let n_rel = causet.causal_matrix().relation_count();
        let n_link = causet.link_matrix().link_count();
        let link_ratio = n_link as f64 / n as f64;

        println!(
            "{:>5}  {:>6}  {:>9}  {:>5}    {:.4}    {:.3}",
            trial + 1, n, n_rel, n_link, f_obs, link_ratio
        );

        ordering_fractions.push(f_obs);
        link_ratios.push(link_ratio);
    }

    let f_mean: f64 = ordering_fractions.iter().sum::<f64>() / n_trials as f64;
    let f_std: f64 = if n_trials > 1 {
        (ordering_fractions
            .iter()
            .map(|&f| (f - f_mean).powi(2))
            .sum::<f64>()
            / (n_trials - 1) as f64)
            .sqrt()
    } else {
        0.0
    };
    let link_mean: f64 = link_ratios.iter().sum::<f64>() / n_trials as f64;

    println!();
    println!("Mean ordering fraction: {:.6} ± {:.6}", f_mean, f_std);
    println!("Mean links per element: {:.3}", link_mean);

    (f_mean, f_std, link_mean)
}

fn print_summary(dim: &str, results: &(f64, f64, f64)) {
    println!(
        "{}    {:.4} ± {:.4}      {:.3}",
        dim, results.0, results.1, results.2
    );
}
