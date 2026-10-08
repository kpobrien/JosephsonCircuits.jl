# Numerical references

These references describe the numerical methods used or adapted in the
package. The implementation notes and tests define the package-specific
choices; a reference to an algorithm does not imply that every detail of
its published implementation is used unchanged.

## Harmonic balance and circuit models

The [home-page bibliography](index.md#References) lists the circuit and
harmonic-balance texts and the published device comparisons.

## Nonlinear and linear solvers

- S. C. Eisenstat and H. F. Walker, “Choosing the Forcing Terms in an
  Inexact Newton Method,” *SIAM Journal on Scientific Computing* 17(1),
  16–32 (1996). [DOI](https://doi.org/10.1137/0917003),
  [authors' technical report](https://softlib.rice.edu/pub/CRPC-TRs/reports/CRPC-TR94463.pdf).
- Y. Saad and M. H. Schultz, “GMRES: A Generalized Minimal Residual
  Algorithm for Solving Nonsymmetric Linear Systems,” *SIAM Journal on
  Scientific and Statistical Computing* 7(3), 856–869 (1986).
  [DOI](https://doi.org/10.1137/0907058).
- D. G. Anderson, “Iterative Procedures for Nonlinear Integral Equations,”
  *Journal of the ACM* 12(4), 547–560 (1965).
  [DOI](https://doi.org/10.1145/321296.321305).
- H. F. Walker and P. Ni, “Anderson Acceleration for Fixed-Point
  Iterations,” *SIAM Journal on Numerical Analysis* 49(4), 1715–1735
  (2011). [DOI](https://doi.org/10.1137/10078356X).
- J. M. Tang, R. Nabben, C. Vuik, and Y. A. Erlangga, “Comparison of
  Two-Level Preconditioners Derived from Deflation, Domain Decomposition
  and Multigrid Methods,” *Journal of Scientific Computing* 39(3), 340–370
  (2009). [DOI](https://doi.org/10.1007/s10915-009-9272-6).

These cover the forcing terms by which Newton–Krylov sets the accuracy of
each linear solve, the GMRES that solves it, the Anderson acceleration of
`QuasiNewton`, and the A-DEF1 form of the optional Floquet deflation.

## Stability

W.-J. Beyn, “An Integral Method for Solving Nonlinear Eigenvalue
Problems,” *Linear Algebra and its Applications* 436(10), 3839–3863
(2012). [DOI](https://doi.org/10.1016/j.laa.2011.03.030).

`ContourIntegral` reduces the moments of the inverse of the perturbation
operator around its circle to the poles inside, as in this paper; see the
[stability guide](stability.md).

## Quantum noise

The noise conventions, symmetrized noise in quanta with the vacuum at half
a photon, and the quantum limit of a phase-preserving amplifier, half a
photon of added noise referred to its input, follow

- C. M. Caves, “Quantum Limits on Noise in Linear Amplifiers,” *Physical
  Review D* 26(8), 1817–1839 (1982).
  [DOI](https://doi.org/10.1103/PhysRevD.26.1817).
- A. A. Clerk, M. H. Devoret, S. M. Girvin, F. Marquardt, and
  R. J. Schoelkopf, “Introduction to Quantum Noise, Measurement, and
  Amplification,” *Reviews of Modern Physics* 82(2), 1155–1208 (2010).
  [DOI](https://doi.org/10.1103/RevModPhys.82.1155),
  [arXiv](https://arxiv.org/abs/0810.4729).

The noise a passive network emits in equilibrium, `(nbar + 1/2)(I - S*S')`,
the noise of a passive block and of `calcCnoise`, is Bosma's theorem:
H. Bosma, *On the Theory of Linear Noisy Systems*, Ph.D. thesis,
Technische Hogeschool Eindhoven (1967).

## Gauss collocation

E. Hairer, *Geometric Numerical Integration*, Lecture 2: Symplectic
integrators (2010), discusses Gauss collocation and preservation of
quadratic invariants. [Lecture notes](https://www.unige.ch/~hairer/poly_geoint/week2.pdf).
These results concern the underlying integration method; algebraic
projection and history interpolation require their own analysis.

## Rational fitting

- B. Gustavsen and A. Semlyen, “Rational Approximation of Frequency Domain
  Responses by Vector Fitting,” *IEEE Transactions on Power Delivery*
  14(3), 1052–1061 (1999).
  [Paper](https://www.sintef.no/globalassets/project/vectfit/vf_paper.pdf).
- B. Gustavsen, “Improving the Pole Relocating Properties of Vector
  Fitting,” *IEEE Transactions on Power Delivery* 21(3), 1587–1592 (2006).
  [Paper](https://www.sintef.no/globalassets/project/vectfit/relaxed_vf_paper.pdf).
- D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
  “Macromodeling of Multiport Systems Using a Fast Implementation of the
  Vector Fitting Method,” *IEEE Microwave and Wireless Components Letters*
  18(6), 383–385 (2008).
  [Paper](https://www.sintef.no/globalassets/project/vectfit/fast_vf_paper.pdf).
- B. Gustavsen, “An Efficient Residue Perturbation Scheme for Passivity
  Enforcement of S-Parameter Rational Models,” *IEEE Transactions on
  Electromagnetic Compatibility* 67(3), 913–920 (2025).
  [Paper](https://doi.org/10.1109/TEMC.2025.3564403).
- C. Eckart and G. Young, “The Approximation of One Matrix by Another of
  Lower Rank,” *Psychometrika* 1(3), 211–218 (1936).
  [Paper](https://doi.org/10.1007/BF02288367).
- J. R. Rice and K. H. Usow, “The Lawson Algorithm and Extensions,”
  *Mathematics of Computation* 22(101), 118–127 (1968).
  [Paper](https://doi.org/10.1090/S0025-5718-1968-0232137-6).

These cover the relocation of common poles, the relaxed normalization,
the per-entry reduction that accelerates a multiport fit, and the
perturbation of the residues that enforces passivity. The last two serve
the order search: a fit of `N` poles is no closer to the samples than
their best approximation of rank `N + 1` (Eckart and Young), so an order
whose bound exceeds the tolerance is not fitted, and Lawson's reweighted
least squares bounds from below the error of every fit with given poles
(Rice and Usow), so an order whose poles cannot meet the tolerance is not
made passive.

## Transfer-matrix norms

S. Boyd and V. Balakrishnan, “A Regularity Result for the Singular Values
of a Transfer Matrix and a Quadratically Convergent Algorithm for Computing
its L-infinity-norm,” *Systems & Control Letters* 15(1), 1–7 (1990).
[Author's page and paper](https://web.stanford.edu/~boyd/papers/sv_of_tf.html).

The package finds the frequencies where a singular value of a rational
realization crosses a level from the Hamiltonian matrix or a pencil, as
in this paper, to validate a realization supplied directly at its
tolerance. See
[`passivityassessment`](@ref JosephsonCircuits.passivityassessment) for interpreting its verdicts.

## Traveling-wave line models

H. W. Dommel, “Digital Computer Solution of Electromagnetic Transients in
Single- and Multiphase Networks,” *IEEE Transactions on Power Apparatus
and Systems* PAS-88(4), 388–399 (1969).
[DOI](https://doi.org/10.1109/TPAS.1969.292459).

This is a historical circuit-transient reference for transmission-line
models based on traveling-wave delay relations. The package's centered
history interpolation, Gauss stages, and discrete adjoint are described
separately in the [transient theory](transienttheory.md#Transmission-lines).
