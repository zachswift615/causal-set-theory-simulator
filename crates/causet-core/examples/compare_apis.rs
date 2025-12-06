//! Compare old vs new Minkowski sprinkling APIs

use causet_core::prelude::*;

fn main() {
    println!("=== Comparing old vs new Minkowski sprinkling ===\n");

    let n_points = 300;
    let seed = 42_u64;

    // Old API: CausalDiamond + original Sprinkler
    let diamond = CausalDiamond::<4>::symmetric(10.0);
    let mut sprinkler = Sprinkler::with_seed(seed, 1.0);
    let old_result = sprinkler.sprinkle_fixed(&diamond, n_points);
    let old_causet = CausalSet::from_sprinkling(old_result);
    let old_action = old_causet.bd_action_4d();

    println!("OLD API (CausalDiamond, τ = 10):");
    println!("  Diamond volume: {:.4}", diamond.volume());
    println!("  Ordering fraction: {:.4}", old_causet.ordering_fraction());
    println!("  BD action: {:.1}", old_action.action);
    println!(
        "  Interval counts: N={}, N1={}, N2={}, N3={}, N4={}",
        old_action.n_elements,
        old_action.interval_counts.n_k(0),
        old_action.interval_counts.n_k(1),
        old_action.interval_counts.n_k(2),
        old_action.interval_counts.n_k(3)
    );

    // New API: Minkowski spacetime with same proper time
    let mink = Minkowski::<4>::causal_diamond(5.0); // half_height = 5 → τ = 10
    let new_result = sprinkle_spacetime(&mink, n_points, seed);
    let new_causet = CausalSet::from_generic_sprinkling(new_result, &mink);
    let new_action = new_causet.bd_action_4d();

    println!("\nNEW API (Minkowski spacetime, half_height = 5, τ = 10):");
    println!("  Minkowski volume: {:.4}", mink.volume());
    println!("  Ordering fraction: {:.4}", new_causet.ordering_fraction());
    println!("  BD action: {:.1}", new_action.action);
    println!(
        "  Interval counts: N={}, N1={}, N2={}, N3={}, N4={}",
        new_action.n_elements,
        new_action.interval_counts.n_k(0),
        new_action.interval_counts.n_k(1),
        new_action.interval_counts.n_k(2),
        new_action.interval_counts.n_k(3)
    );

    println!("\n=== Dimension recovery check ===");
    println!("Old API estimated dimension: {:.2}", old_causet.estimated_dimension().dimension);
    println!("New API estimated dimension: {:.2}", new_causet.estimated_dimension().dimension);

    // Verify all points are inside the diamond
    println!("\n=== Diamond containment check ===");
    let half_h = 5.0_f64;

    let old_outside: Vec<_> = old_causet.points().iter()
        .filter(|p| {
            let t = p.coords[0].abs();
            let r: f64 = (1..4).map(|i| p.coords[i].powi(2)).sum::<f64>().sqrt();
            t + r > half_h + 0.01  // small tolerance
        })
        .collect();

    let new_outside: Vec<_> = new_causet.points().iter()
        .filter(|p| {
            let t = p.coords[0].abs();
            let r: f64 = (1..4).map(|i| p.coords[i].powi(2)).sum::<f64>().sqrt();
            t + r > half_h + 0.01
        })
        .collect();

    println!("Old API: {} points outside diamond", old_outside.len());
    println!("New API: {} points outside diamond", new_outside.len());

    // Print a few if any are outside
    for p in new_outside.iter().take(3) {
        let t = p.coords[0].abs();
        let r: f64 = (1..4).map(|i| p.coords[i].powi(2)).sum::<f64>().sqrt();
        println!("  Outside: t={:.2}, r={:.2}, sum={:.2} (should be <= {:.2})", t, r, t+r, half_h);
    }

    // Check time distribution
    println!("\n=== Time distribution ===");
    let old_mean_t: f64 = old_causet.points().iter().map(|p| p.coords[0].abs()).sum::<f64>() / n_points as f64;
    let new_mean_t: f64 = new_causet.points().iter().map(|p| p.coords[0].abs()).sum::<f64>() / n_points as f64;
    println!("Old API mean |t|: {:.3}", old_mean_t);
    println!("New API mean |t|: {:.3}", new_mean_t);

    // Check spatial spread
    let old_mean_r: f64 = old_causet.points().iter().map(|p| {
        (1..4).map(|i| p.coords[i].powi(2)).sum::<f64>().sqrt()
    }).sum::<f64>() / n_points as f64;
    let new_mean_r: f64 = new_causet.points().iter().map(|p| {
        (1..4).map(|i| p.coords[i].powi(2)).sum::<f64>().sqrt()
    }).sum::<f64>() / n_points as f64;
    println!("Old API mean r: {:.3}", old_mean_r);
    println!("New API mean r: {:.3}", new_mean_r);
}
