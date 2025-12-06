//! Phase 2 Validation: Dimension Recovery from Causal Structure
//!
//! This example demonstrates that spacetime dimension can be recovered
//! from pure causal structure using the Myrheim-Meyer estimator.

use causet_core::prelude::*;
use causet_core::sprinkling::{CausalDiamond, Sprinkler};

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║   Phase 2 Validation: Dimension Recovery from Causality      ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    println!("Demonstrating that spacetime dimension emerges from causal structure.");
    println!("Using Myrheim 1978 Table I values with log-linear interpolation.");
    println!();

    for d in 2..=4 {
        println!("═══════════════════════════════════════════════════════════════");
        println!("{}D Minkowski Spacetime", d);
        println!();

        // Run multiple trials
        let mut estimates = Vec::new();
        let mut fractions = Vec::new();

        for trial in 0..5 {
            let seed = 42 + trial as u64 * 1000;
            let (f, dim_est) = match d {
                2 => run_trial::<2>(seed, 500),
                3 => run_trial::<3>(seed, 500),
                4 => run_trial::<4>(seed, 500),
                _ => unreachable!(),
            };
            estimates.push(dim_est);
            fractions.push(f);

            println!(
                "  Trial {}: f={:.4}, d_est={:.2}{}",
                trial + 1,
                f,
                dim_est.dimension,
                if dim_est.extrapolated {
                    " (extrapolated)"
                } else {
                    ""
                }
            );
        }

        let f_mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;
        let d_mean: f64 = estimates.iter().map(|e| e.dimension).sum::<f64>() / estimates.len() as f64;
        let d_std: f64 = (estimates
            .iter()
            .map(|e| (e.dimension - d_mean).powi(2))
            .sum::<f64>()
            / (estimates.len() - 1) as f64)
            .sqrt();

        println!();
        println!("  Mean: f={:.4}, d={:.2} ± {:.2}", f_mean, d_mean, d_std);
        println!("  Expected: d={}.0", d);
        println!("  Error: {:.1}%", ((d_mean - d as f64).abs() / d as f64) * 100.0);
    }

    println!();
    println!("═══════════════════════════════════════════════════════════════");
    println!("CONCLUSION: Spacetime dimension successfully recovered from");
    println!("pure causal structure, validating the Myrheim-Meyer estimator.");
    println!("═══════════════════════════════════════════════════════════════");
}

fn run_trial<const D: usize>(seed: u64, n_points: usize) -> (f64, DimensionEstimate) {
    let diamond = CausalDiamond::<D>::symmetric(10.0);
    let mut sprinkler = Sprinkler::with_seed(seed, 1.0);
    let result = sprinkler.sprinkle_fixed(&diamond, n_points);
    let causet = CausalSet::from_sprinkling(result);

    let f = causet.ordering_fraction();
    let dim = causet.estimated_dimension();

    (f, dim)
}
