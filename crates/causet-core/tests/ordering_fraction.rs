//! Integration tests for ordering fraction computation.

use causet_core::prelude::*;

#[test]
fn test_chain_has_unit_ordering_fraction() {
    // A chain: all pairs causally related
    let points: Vec<Point2D> = (0..10)
        .map(|i| Point2D::new_2d(i as f64, 0.0, i))
        .collect();

    let causet = CausalSet::from_points(points);
    assert!((causet.ordering_fraction() - 1.0).abs() < 0.001);
}

#[test]
fn test_antichain_has_zero_ordering_fraction() {
    // An antichain: no pairs causally related (all at same time)
    let points: Vec<Point2D> = (0..10)
        .map(|i| Point2D::new_2d(0.0, i as f64, i))
        .collect();

    let causet = CausalSet::from_points(points);
    assert!((causet.ordering_fraction() - 0.0).abs() < 0.001);
}

#[test]
fn test_known_diamond_geometry() {
    // Diamond: bottom → two middle → top
    // 5 relations out of 6 pairs = f ≈ 0.833
    let points = vec![
        Point2D::new_2d(0.0, 0.0, 0),
        Point2D::new_2d(1.0, -0.5, 1),  // Timelike from (0,0)
        Point2D::new_2d(1.0, 0.5, 2),   // Timelike from (0,0)
        Point2D::new_2d(2.0, 0.0, 3),   // Timelike from all
    ];

    let causet = CausalSet::from_points(points);
    let stats = causet.statistics();

    // Verify: (0,0)→(1,-0.5), (0,0)→(1,0.5), (0,0)→(2,0), (1,-0.5)→(2,0), (1,0.5)→(2,0)
    // = 5 relations
    assert_eq!(stats.n_relations, 5);
    assert!((causet.ordering_fraction() - 5.0/6.0).abs() < 0.01);
}

#[test]
fn test_sprinkling_produces_intermediate_fraction() {
    // Random sprinkling should give 0 < f < 1
    let diamond = CausalDiamond::<2>::symmetric(4.0);
    let mut sprinkler = Sprinkler::with_seed(12345, 50.0);
    let result = sprinkler.sprinkle_diamond(&diamond);
    let causet = CausalSet::from_sprinkling(result);

    let f = causet.ordering_fraction();
    assert!(f > 0.1, "Ordering fraction too low: {}", f);
    assert!(f < 0.99, "Ordering fraction too high: {}", f);
}
