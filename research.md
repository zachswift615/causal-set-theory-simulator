# Building causal set quantum gravity simulations from scratch

Causal set theory offers a discrete approach to quantum gravity where spacetime emerges from a partially ordered set of events, with the causal ordering encoding geometric information. This guide provides the complete technical specifications needed to implement working simulations, covering all algorithms, formulas, and data structures.

## Sprinkling points into Lorentzian manifolds

The foundation of causal set simulations is **Poisson sprinkling**—placing points randomly into spacetime according to a Poisson process that preserves Lorentz invariance (unlike regular lattices). The probability of placing exactly k points in a region R with spacetime volume V_R follows:

```
P(N(R) = k) = (ρV_R)^k · e^(-ρV_R) / k!
```

The **density parameter ρ** has units of inverse spacetime volume and relates to the discreteness scale ℓ via **ρ = ℓ^(-d)** in d dimensions. Physically, ℓ is presumed to be the Planck length (≈10⁻³⁵ m), meaning each causal set element corresponds to roughly one Planck volume. For simulations, typical values range from ρ = 1 to ρ = 1000 in natural units.

**Minkowski spacetime sprinkling** uses the metric ds² = -dt² + dx² (in 2D) with volume element dV = dt dx. For a causal diamond (Alexandrov interval) between points p and q, the algorithm samples uniformly from a bounding box and rejects points outside the diamond. The diamond volume in d dimensions is **V_d = C_d · τ^d** where τ is proper time separation, with coefficients **C₂ = 1/2, C₃ = π/12, and C₄ = π/24**. The general formula is C_d = π^((d-1)/2) / [2^(d-1) · d · Γ((d+1)/2)].

For **curved spacetimes like de Sitter**, the volume element becomes position-dependent. In global coordinates with metric ds² = -ℓ² dt² + ℓ² cosh²(t) dθ², importance sampling handles the time-varying volume element through inverse CDF methods. The **Schwarzschild spacetime** requires Eddington-Finkelstein Original (EFO) coordinates where constant-t* hypersurfaces remain spacelike everywhere, enabling clean volume integration via dV = r² sin(θ) dt* dr dθ dφ.

## Computing causal relations efficiently

Given two points with coordinates, the causal relation test in Minkowski spacetime checks whether **Δs² = -(t_y - t_x)² + |x_y - x_x|² ≤ 0 AND t_y > t_x**. For Schwarzschild spacetime, the algorithm distinguishes exterior and interior regions:

**Outside horizon (r₁, r₂ > 2M):** For ingoing null geodesics, check dt ≥ r₁ - r₂; for outgoing, check dt ≥ r₂ - r₁ + 4M·ln[(r₂-2M)/(r₁-2M)].

**Inside horizon (r₁, r₂ < 2M):** The constraint becomes r₁ - r₂ ≤ dt ≤ r₂ - r₁ + 4M·ln[(2M-r₂)/(2M-r₁)].

Computing the full **N×N causal matrix** requires O(N²) comparisons minimum. Key optimizations include: sorting points by time coordinate first (only checking pairs where t_j > t_i), spatial pruning with bounding boxes, and GPU parallelization since each pair check is independent.

**Data structures** trade off between space and query efficiency. The full adjacency matrix C[i,j] = 1 if p_i ≺ p_j uses O(N²) bits but offers O(1) lookup. The **link matrix** stores only direct causal connections (transitive reduction)—far sparser for manifold-like causets. DAG adjacency lists using parents/children arrays achieve O(N + E) space where E is link count. The **FastBitset** implementation uses compressed-bit storage with SIMD instructions for orders-of-magnitude speedup.

For transitive reduction (computing links from the full causal matrix), Hsu's O(N³) algorithm iteratively removes redundant edges, while the Aho-Garey-Ullman approach uses matrix multiplication achieving O(N^2.373) complexity.

## The Myrheim-Meyer dimension estimator

This elegant technique extracts spacetime dimension from causal structure alone. For N elements in an order interval with R causally related pairs, the **ordering fraction** is:

```
f = 2R / [N(N-1)]
```

The fundamental relationship for d-dimensional Minkowski spacetime:

```
f(d) = Γ(d+1) · Γ(d/2) / [4 · Γ(3d/2)]
```

**IMPORTANT:** In this formula, **d = number of SPATIAL dimensions**, not spacetime dimensions!

From Myrheim's original 1978 CERN preprint (TH-2538), Table I:

| Spacetime | Spatial d | Ordering Fraction |
|-----------|-----------|-------------------|
| 1+1D      | 1         | 1/2 = 0.500       |
| 2+1D      | 2         | 8/35 ≈ 0.229      |
| 3+1D      | 3         | 1/10 = 0.100      |

**Note:** Some secondary sources incorrectly cite f(3D)≈0.424 and f(4D)≈0.333 — these values are wrong and do not appear in the primary literature. Always verify against Myrheim 1978.

The function monotonically decreases, so inverting numerically gives a unique dimension estimate from observed ordering fraction.

The **midpoint-scaling variant** finds the element z maximizing min(|I(x,z)|, |I(z,y)|), then computes **d ≈ ln(N/N_mid) / ln(2)**. This approach is invariant under coarse-graining and computationally simpler. Statistical uncertainty scales as **δd/d ≈ 1/√N**, with systematic bias toward lower dimensions for N < 100.

```python
def myrheim_meyer_dimension(causal_matrix):
    N = len(causal_matrix)
    R = np.sum(causal_matrix)  # Count all related pairs
    f_obs = 2.0 * R / (N * (N - 1))
    
    def f_of_d(d):
        return gamma(d+1) * gamma(d/2) / (4 * gamma(3*d/2)) - f_obs
    
    return brentq(f_of_d, 1.5, 10.0)  # Numerical inversion
```

## Benincasa-Dowker action with exact coefficients

The BD action defines discrete scalar curvature through a causal set d'Alembertian. The complete formula for d dimensions:

```
(1/ℏ) S^(d)_BDG(C) = ζ_d · [N + (β_d/α_d) · Σ_{i=1}^{n_d} C^(d)_i · N_i]
```

where N_i counts order intervals of cardinality (i+1), and the coefficients are:

| Dimension | α_d | β_d | n_d | C_k coefficients |
|-----------|-----|-----|-----|------------------|
| 2 | -2 | 4 | 3 | (1, -2, 1) |
| 3 | -3π/√(2π) | 6π/√(2π) | 3 | (1, -3, 2) |
| 4 | -2√6 | 4√6 | 4 | (1, -9, 16, -8) |

The **4D action** simplifies to: **S^(4) ∝ N - N_1 + 9N_2 - 16N_3 + 8N_4** (with ℓ = ℓ_p). Here N_k counts pairs (x,y) where x ≺ y and exactly k elements lie strictly between them. The continuum limit recovers the Einstein-Hilbert action plus a boundary term: **lim_{ρ→∞} ⟨S_ρ(M)⟩ = (1/ℓ_p^{d-2})[½∫_M d^dx √(-g)R + Vol_{d-2}(J)]**.

Computing the action requires counting intervals of each size—an O(N³) operation with the best known classical algorithm, reducible to O(N^2.8) with sufficient memory via matrix multiplication.

## Classical Sequential Growth dynamics

CSG provides a stochastic model for growing causal sets element-by-element. The **transitive percolation** model uses a single parameter p ∈ [0,1]:

1. Start with empty set or single element
2. Add new element y to existing n elements
3. For each element x_i: with probability p, create link x_i ≺ y
4. Apply transitivity: if x_i ≺ x_j and x_j ≺ y, then x_i ≺ y
5. Repeat

The transition probability follows **P(C_n → C_{n+1}) = p^m · (1-p)^{n-k}** where m is the number of maximal elements in the precursor set and k is the precursor set size. The general CSG model satisfies two constraints: **discrete general covariance** (path independence—probability of reaching a final causet is independent of labeling order) and **Bell causality** (birth affected only by causal past, not spacelike elements).

The coupling constants {t₀, t₁, t₂, ...} parameterize allowed models through q_n = Σ_{j=0}^n t_j · (n choose j). Physical interpretation: small p produces "thin" causal sets with few relations; large p produces "thick" ones approaching total orders. **Posts** (elements remaining maximal and related to all future elements) divide causal evolution into cosmological epochs.

## Observables for measuring geometry

**Proper time estimation** uses the longest chain between causally related points: τ(x,y) ∝ max{chain lengths from x to y}. Calibrated against sprinkling density, **τ = C_d · L / ρ^{1/d}** where L is longest chain length.

**Spatial distance** between spacelike-separated events uses the **causal overlap method**: for events a, b with common future event c, compute the overlap of their causal diamonds normalized by relevant volumes. This achieves Planck-scale accuracy in Minkowski spacetimes.

**Manifold-likeness tests** distinguish physical causal sets from random partial orders. The Kleitman-Rothschild theorem shows asymptotically almost all partial orders are 3-layer structures (definitely not manifold-like). Tests include: consistent dimension estimation via Myrheim-Meyer, BD action suppression of pathological orders in path integrals, and stable homology under thickened antichain constructions.

**Thickened antichains** provide spatial hypersurface identification: starting from an inextendible antichain A, the future-volume thickening A_v = {x : |past(x) ∩ future(A)| ≤ v} creates convex regions carrying topological information recoverable via nerve complex homology.

## Available implementations and repositories

**Python-causets** (github.com/c-minz/Python-causets) offers the most accessible starting point: Python 3.8+ with mypy typing, supporting sprinkling on flat/de Sitter/black hole spacetimes, plotting with light-cones, and animations. BSD 3-Clause licensed.

**CausalSetTools** (Computer Physics Communications) provides high-performance C/C++ with CUDA GPU acceleration, achieving orders of magnitude speedup through FastBitset structures and AVX2 optimizations. Essential for N > 10,000 elements.

**energetic-causal-set-simulation** (github.com/gfmio) implements the Cortês-Smolin energetic causal set model in Python, useful for understanding event-based dynamics.

David Rideout's UCSD resources (mathweb.ucsd.edu/~drideout) include poscau diagrams up to 5 elements and links to foundational theses.

| Causal Set Size | Feasible Operations | Typical Time |
|-----------------|---------------------|--------------|
| N < 1,000 | Full analysis, all estimators | Minutes |
| N ~ 10,000 | Sprinkling, basic topology | Hours |
| N ~ 100,000 | Sprinkling with GPU | Hours-Days |
| N > 1,000,000 | Generation only | Requires HPC |

## Essential references for implementation

The **foundational papers** are Bombelli-Lee-Meyer-Sorkin (1987) "Space-time as a causal set" (Phys. Rev. Lett. 59, 521) defining axioms and sprinkling, and Myrheim's 1978 CERN preprint TH-2538 introducing statistical geometry. Meyer's 1988 MIT PhD thesis (hdl.handle.net/1721.1/14328) contains the complete dimension estimator derivation.

For **dynamics**, Rideout-Sorkin (2000) "A Classical Sequential Growth Dynamics" (Phys. Rev. D 61, 024002; arXiv:gr-qc/9904062) provides algorithmic CSG descriptions. The **Benincasa-Dowker action** appears in Phys. Rev. Lett. 104, 181301 (2010; arXiv:1001.2725), with arbitrary-dimension extensions in Dowker-Glaser (2013; arXiv:1305.2588) and closed-form expressions in Glaser (2014; arXiv:1311.1701).

**Surya's Living Reviews in Relativity article** (2019; arXiv:1903.11544) is the essential comprehensive reference, covering all major topics semi-pedagogically with extensive bibliography. For topology recovery, Major-Rideout-Surya (2007) "On recovering continuum topology" (J. Math. Phys. 48, 032501) details thickened antichain homology.

## Conclusion

Building a working causal set simulation requires implementing: (1) Poisson sprinkling with correct volume measures for target spacetimes, (2) efficient O(N²) causal relation computation with sparse storage, (3) dimension estimation via Myrheim-Meyer inversion, (4) BD action through interval counting with coefficients matched to dimension, and (5) optionally CSG dynamics via transitive percolation. The Python-causets package provides an excellent starting template, while CausalSetTools enables scaling to physically interesting sizes. The critical insight enabling all geometric recovery is that causal structure plus cardinality approximately determines spacetime geometry—the core conjecture underlying the entire program.