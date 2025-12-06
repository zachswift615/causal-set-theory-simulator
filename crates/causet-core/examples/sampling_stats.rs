//! Statistical comparison of old vs new sampling

use causet_core::prelude::*;

fn main() {
    println!("=== Statistical sampling comparison (20 trials each) ===\n");

    let n_points = 500;
    let n_trials = 20;

    let mut old_ord_fracs = Vec::new();
    let mut new_ord_fracs = Vec::new();
    let mut old_actions = Vec::new();
    let mut new_actions = Vec::new();

    for seed in 0..n_trials {
        // Old API
        let diamond = CausalDiamond::<4>::symmetric(10.0);
        let mut sprinkler = Sprinkler::with_seed(seed as u64, 1.0);
        let old_result = sprinkler.sprinkle_fixed(&diamond, n_points);
        let old_causet = CausalSet::from_sprinkling(old_result);
        old_ord_fracs.push(old_causet.ordering_fraction());
        old_actions.push(old_causet.bd_action_4d().action);

        // New API
        let mink = Minkowski::<4>::causal_diamond(5.0);
        let new_result = sprinkle_spacetime(&mink, n_points, seed as u64);
        let new_causet = CausalSet::from_generic_sprinkling(new_result, &mink);
        new_ord_fracs.push(new_causet.ordering_fraction());
        new_actions.push(new_causet.bd_action_4d().action);
    }

    let old_mean_f: f64 = old_ord_fracs.iter().sum::<f64>() / n_trials as f64;
    let new_mean_f: f64 = new_ord_fracs.iter().sum::<f64>() / n_trials as f64;
    let old_mean_s: f64 = old_actions.iter().sum::<f64>() / n_trials as f64;
    let new_mean_s: f64 = new_actions.iter().sum::<f64>() / n_trials as f64;

    let old_std_f: f64 = (old_ord_fracs.iter().map(|x| (x - old_mean_f).powi(2)).sum::<f64>() / (n_trials - 1) as f64).sqrt();
    let new_std_f: f64 = (new_ord_fracs.iter().map(|x| (x - new_mean_f).powi(2)).sum::<f64>() / (n_trials - 1) as f64).sqrt();
    let old_std_s: f64 = (old_actions.iter().map(|x| (x - old_mean_s).powi(2)).sum::<f64>() / (n_trials - 1) as f64).sqrt();
    let new_std_s: f64 = (new_actions.iter().map(|x| (x - new_mean_s).powi(2)).sum::<f64>() / (n_trials - 1) as f64).sqrt();

    println!("Ordering fraction (theoretical for 4D: 0.100):");
    println!("  Old API: {:.4} ± {:.4}", old_mean_f, old_std_f);
    println!("  New API: {:.4} ± {:.4}", new_mean_f, new_std_f);

    println!("\nBD action (expected ~0 for flat spacetime):");
    println!("  Old API: {:.1} ± {:.1}", old_mean_s, old_std_s);
    println!("  New API: {:.1} ± {:.1}", new_mean_s, new_std_s);

    println!("\nExpected fluctuation (√N): {:.1}", (n_points as f64).sqrt());
}
