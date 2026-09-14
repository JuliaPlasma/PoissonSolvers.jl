# Release Notes

All notable changes to PoissonSolvers.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31. 11 versions were released before it, the most recent
`v0.3.7`, and none of them are written up here: the record of that history is `git log` and
the tags. That gap is deliberate, and is named rather than
reconstructed, because a changelog assembled after the fact loses exactly the reasoning that
makes it worth keeping. The record proper begins with the section below.

## [Unreleased] — targeting 0.4.0

### New Features

- **`DirichletBasisSpline(domain, order, ncells)` is a new backend** for non-singular solves of
  ``-\phi'' = \rho`` with homogeneous Dirichlet boundaries. The basis is recombined so every basis
  function vanishes at both ends, and the boundary condition holds exactly rather than
  approximately. Measured convergence is at the order of the basis (rate 3.0 for order 3, 5.0 for
  order 5), same as the periodic solver.

- **Aqua.jl is now part of the test suite**, guarding against type piracy, undefined exports, stale
  dependencies and the other faults a behaviour-only suite cannot see. The suite grew from about
  20 assertions to 161, adding `update!`, `rhs`, in-place `solve!`, derivative evaluation, the
  zero-source `Potential` constructor, convergence rates, order-vs-degree checks, error paths, and
  type-stability and allocation gates.

### Bug Fixes

- **`PeriodicBasisFFT` is now defined and exported**. The function existed but was named
  `PeriodicBasisFFTW`, leaving the exported name undefined — the fault Aqua's `undefined_exports`
  now guards against.

- **Derivative evaluation no longer rebuilds the basis on every call.** It was reconstructing the
  entire spline basis and coefficient vector for every evaluation (from BSplineKit's `Derivative(1)
  * spline`). Now it is a stateless lookup, which is visible to downstream code because the
  evaluation happens once per particle per stage.

- **`solve!` is now allocation-free on all backends.** The spline solver was rebuilding an FFT
  factorization on every call, and the grid solver was rebuilding FFT plans and three temporary
  arrays. Transforms and factorizations are now constructed once when the solver is built.
  Measured: `solve!` and `update!` allocate 0 bytes on the periodic spline, Dirichlet spline and
  grid backends alike.

- **The grid solver's constant mode no longer divides by zero.** The k = 0 Fourier coefficient was
  computed by dividing by zero and then overwritten; the symbol now carries an exact zero there.

- **`FastGaussQuadrature` was a declared dependency that no code used**, and is removed. Aqua's
  `stale_deps` did not catch it and could not have: it was also a transitive dependency of
  BSplineKit, so it was present in the loaded closure and looked used. Dropping BSplineKit is what
  makes it visible.

### Breaking Changes

- **The B-spline backend now uses SimpleSplines.jl instead of BSplineKit.jl**, and the associated
  dependencies are gone from `[deps]`: BSplineKit, ToeplitzMatrices, SparseArrays and
  FastGaussQuadrature. SimpleSplines 0.1 is added. The Julia floor remains 1.10.

- **Spline solvers are renamed** to describe the method rather than the library: `PeriodicBasisBSplineKit`
  → `PeriodicBasisSpline`, and `PoissonSolverBSplineKit` → `PoissonSolverSpline`. Grid solvers are
  unchanged.

- **The derivative API matches SimpleSplines' spelling.** `p(x, Derivative(1))` becomes `p(x, 1)`,
  and `derivative(p, d)` returns a callable for broadcasting. The re-exported `BSplineKit.Derivative`
  is gone. PoissonSolvers now re-exports `basis`, `coefficients` and `derivative` from SimpleSplines,
  so a caller with `using PoissonSolvers` reaches them the same way they would reach a
  `SimpleSplines.Spline`.

- **`Potential` now has three fields** (`potential`, `solver`, `rhs`) instead of five. Access the
  basis and coefficients through the functions `basis(p)` and `coefficients(p)` rather than the
  struct fields `p.basis` and `p.coefficients`, matching SimpleSplines' idiom. Both were already
  stored inside the solution, so this is renaming access rather than adding storage.

- **`evalsolution` is removed.** Build a `Spline(basis, coeffs)` from the components or use the
  `Potential` functor.

- **`PoissonSolverSpline` no longer accepts arbitrary non-periodic bases** through a branch that
  left two struct fields undefined and could not have worked. It takes a periodic or a Dirichlet
  basis only.

- **Minimum Julia is now 1.10**, raised from the declared 1.8. 1.10 is the LTS and the floor
  across the whole tree; 1.8 was declared but never tested and would not resolve against the
  current dependency versions. CI now derives its lower matrix entry from this field, so a
  declared floor that nobody tests is no longer possible.

## Open Issues

- **Evaluating a spline potential at a scalar point allocates 96 bytes.** This is not in this
  package: `SimpleSplines.evaluate(basis, coefficients, x)` takes one `local_width` buffer per
  scalar call, because the degree is a struct field rather than a type parameter. The fix belongs
  upstream: `SimpleSplines.evaluate_all!`, which takes a caller-supplied buffer, allocates nothing.
  See `scripts/measure_allocations.jl` for both figures side by side. The grid backend evaluates
  with no allocation.
