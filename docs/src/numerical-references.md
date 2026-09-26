# Numerical references

These references describe the numerical methods used or adapted in the
package. The implementation notes and tests define the package-specific
choices; a reference to an algorithm does not imply that every detail of
its published implementation is used unchanged.

## Harmonic balance and circuit models

The [home-page bibliography](index.md#References) lists the circuit and
harmonic-balance texts and the published device comparisons.

For adaptive linear-solve accuracy in Newton–Krylov, see S. C. Eisenstat
and H. F. Walker, “Choosing the Forcing Terms in an Inexact Newton Method,”
*SIAM Journal on Scientific Computing* 17(1), 16–32 (1996).
[DOI](https://doi.org/10.1137/0917003),
[authors' technical report](https://softlib.rice.edu/pub/CRPC-TRs/reports/CRPC-TR94463.pdf).

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

These cover common-pole relocation, relaxed normalization, and the
per-entry reduction used to accelerate multiport fitting.

## Transfer-matrix norms

S. Boyd and V. Balakrishnan, “A Regularity Result for the Singular Values
of a Transfer Matrix and a Quadratically Convergent Algorithm for Computing
its L-infinity-norm,” *Systems & Control Letters* 15(1), 1–7 (1990).
[Author's page and paper](https://web.stanford.edu/~boyd/papers/sv_of_tf.html).

The package uses a level-set/pencil norm calculation for rational-model
validation, with explicit numerical tolerances. See
[`passivityassessment`](@ref JosephsonCircuits.passivityassessment) for interpreting its reported bounds.

## Traveling-wave line models

H. W. Dommel, “Digital Computer Solution of Electromagnetic Transients in
Single- and Multiphase Networks,” *IEEE Transactions on Power Apparatus
and Systems* PAS-88(4), 388–399 (1969).
[DOI](https://doi.org/10.1109/TPAS.1969.292459).

This is a historical circuit-transient reference for transmission-line
models based on traveling-wave delay relations. The package's centered
history interpolation, Gauss stages, and discrete adjoint are described
separately in the [transient theory](transienttheory.md#Transmission-lines).
