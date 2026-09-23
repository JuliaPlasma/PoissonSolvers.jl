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

## [Unreleased] — targeting 0.6.0

### New Features

- **Matrix-free solver backend `PoissonSolverMatrixFree`, selected by
  `FiniteDifferenceBasis(domain, ngrid; order = 2 | 4)`.** The solver solves
  `(-Δ + R)ϕ = (1 - R)ρ` by conjugate gradients on a uniform periodic grid,
  where `R` projects onto the constant nullspace and the periodic Laplacian is
  the central difference stencil `_apply_Δₓ!` (second order) or
  `_apply_Δₓ₄!` (fourth order). The solver holds iteration buffers as fields,
  so `solve!` from a vector allocates zero bytes; given a function, it allocates
  the vector it samples first. Returns the mean-free solution, the same as
  `PoissonSolverFFT`. The iteration stops once the residual norm falls to
  `reltol` times that of `(1 - R)ρ`, and throws `ErrorException` if `maxiter`
  iterations do not reach it; `solve!` throws `DimensionMismatch` on a vector of
  the wrong length. `FiniteDifferenceBasis` and `PoissonSolverMatrixFree` are
  exported.

- **The stencils `_apply_Δₓ!`, `_apply_Δₓ₄!` and `_apply_Rₓ!` are exported**
  as the operators of the matrix-free backend. They write into a caller-supplied
  vector and build no matrix. They arrive from
  `ReducedBasisMethods/src/gridbased/poisson.jl`; the velocity-moment stencils
  that shared that file went to `VlasovMethods`, which is where a `∫dv` belongs.
  `_apply_Δₓ₄!` and `_apply_Lₓ₄!` differ from the ReducedBasisMethods copies:
  those used (5, −32, 54, −32, 5)/12h², which is second order (measured error on
  sin(2πx): 1.97, 0.503, 0.127, 0.0317 at n = 16, 32, 64, 128); PoissonSolvers
  uses the fourth-order (−1, 16, −30, 16, −1)/12h². `_apply_Lₓ₄!`, the
  regularised `-Δ + R`, is not exported; the backend composes the exported
  stencils instead. Reach it as `PoissonSolvers._apply_Lₓ₄!` if needed.

## [0.5.0] — 2026-09-14

### Breaking Changes

- **`PoissonSolverSpline` no longer regularises the periodic stiffness matrix by a rank-one
  shift.** Constants lie in the kernel of the periodic stiffness matrix, so it is singular. The
  solver now asks `SimpleSplines` for the deflated solve directly, with
  `mass_operator(stiffness_matrix(q), b; kernel = :project)`, which needs `SimpleSplines 0.2`.
  The answers are unchanged: the deflated solve reproduces the shifted one to round-off and
  still returns the mean-free solution. What changes is the cost and the condition number.

  The rank-one shift was a scalar added to every entry, which made a banded assembly
  structurally full. At degree 4 on 2048 cells the shifted matrix alone measured
  **67 125 416 B**, and a whole `PoissonSolverSpline` **72 551 536 B**. The same solver now
  measures **5 426 160 B**, thirteen times smaller, and what it holds is the sparse assembly:
  18 432 stored entries rather than 4 194 304. Construction was O(N²) in both time and storage
  where the assembly itself is O(N).

  The shift also chose a scale it had no basis for: `inv(N)` is an absolute constant while the
  spectrum of the stiffness matrix scales with the mesh. On a domain of length 2π·10⁻³ that
  raised the condition number of the shifted matrix by a factor of 98 over the mean-free
  spectrum. The deflation works on the mean-free subspace itself and has no such scale. This
  is the same choice `PoissonSolverFFT` already makes, where the `k = 0` factor is zero rather
  than `1/k²`.

  The claims above are established by `scripts/verify_kernel_projection.jl`, which is new on
  this branch and checks the deflation against the shift on both representations.

  This resolves the *Open Issues* entry **A periodic spline solver stores an ``n \times n``
  matrix**, filed upstream as
  [SimpleSplines.jl#10](https://github.com/JuliaDEC/SimpleSplines.jl/issues/10). The entry is
  dropped from that section below.

- **`PoissonSolvers.regularise` and `PoissonSolvers.meanfree!` are removed.** Both were internal
  and unexported.

- **`solve!` checks the lengths it was given**, and says how many degrees of freedom the solver
  has when they disagree. The shift's broadcast used to raise that `DimensionMismatch` as a side
  effect, so one corner changes: a result vector of the wrong length, with a right-hand side of
  the right one, used to reach the planned transforms as `ArgumentError: FFTW plan applied to
  wrong-size output`. It is a `DimensionMismatch` now, like every other mismatch.

- **`SimpleSplines` moves from `"0.1"` to `"0.2"`**, for the `kernel` keyword.

## [0.4.0] — 2026-09-14

### New Features

- **`DirichletBasisSpline(domain, order, ncells)` is a new backend** for non-singular solves of
  ``-\phi'' = \rho`` with homogeneous Dirichlet boundaries. The basis is recombined so every basis
  function vanishes at both ends, and the boundary condition holds exactly rather than
  approximately. Measured convergence is at the order of the basis (rate 3.0 for order 3, 5.0 for
  order 5), same as the periodic solver.

- **Aqua.jl is now part of the test suite**, guarding against type piracy, undefined exports, stale
  dependencies and the other faults a behaviour-only suite cannot see. The suite grew from about
  20 assertions to 173, adding `update!`, `rhs`, in-place `solve!`, derivative evaluation, the
  zero-source `Potential` constructor, convergence rates, order-vs-degree checks, error paths, and
  type-stability and allocation gates.

### Bug Fixes

- **The undefined export `PeriodicBasisFFT` is removed.** The name was exported but never defined,
  so every use of it raised `UndefVarError` — the fault Aqua's `undefined_exports` now guards
  against. The defined name was `PeriodicBasisFFTW`, an unexported alias, and it is removed with
  it. `FFTWBasis` is the grid basis constructor, exported and documented.

- **Derivative evaluation no longer rebuilds the basis on every call.** It was reconstructing the
  entire spline basis and coefficient vector for every evaluation (from BSplineKit's `Derivative(1)
  * spline`). Now it is a stateless lookup, which is visible to downstream code because the
  evaluation happens once per particle per stage.

- **`solve!` from a vector is now allocation-free on all backends.** The spline solver was
  rebuilding an FFT factorization on every call, and the grid solver was rebuilding FFT plans and
  three temporary arrays. Transforms and factorizations are now constructed once when the solver
  is built. Measured: `solve!` and `update!` allocate 0 bytes on the periodic spline, Dirichlet
  spline and grid backends alike. Solving with a *function* still allocates the right-hand side it
  samples first — 2080 B on the periodic spline, 2832 B on the Dirichlet one and 704 B on the grid
  at the sizes measured. Fill `rhs(p)` and call `update!(p)` where that matters.

  The scratch that buys this lives in the solver, so **a solver is not reentrant**: two tasks must
  not call `solve!` on one solver, even with distinct result vectors. Give each task its own
  solver. Both solver docstrings say so. Previously every call allocated its own scratch, so this
  is a new constraint rather than one that was always there unstated.

- **`Potential(basis)` no longer forces `Float64`.** The zero right-hand side of the
  single-argument constructor was built as `zeros(ndofs(b))` whatever the basis was, and the two
  backends failed differently. On a `Float32` grid basis it reached the transform as a
  `Vector{Float64}` and threw a `MethodError` from `mul!`. On a `Float32` spline basis it
  constructed without complaint and computed the whole potential at `Float64` — the quieter fault
  and the worse one, since nothing said the requested precision had been discarded. The right-hand
  side now follows `eltype(basis)`, and `Base.eltype` is defined on `FFTWBasis` — it fell back to
  `Any` before, while the SimpleSplines bases already answered it.

- **The regularising shift keeps the element type of the stiffness matrix.** Written as
  `inv(size(S, 1))` the shift was a `Float64` scalar, so it promoted a `Float32` matrix and with it
  the mass operator, which then missed the `MassOperator{DT}` bound on the solver field. No caller
  can reach this today, because a `Float32` periodic basis fails earlier — see `## Open Issues` —
  but the arithmetic is now right whatever the precision.

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

- **`Potential`'s first type parameter is now the solution type, not the basis type.** A method
  dispatching on `Potential{<:SomeBasis}` no longer matches, and because it simply stops being
  called rather than erroring, it fails silently. Dispatch on the solution or on `Potential`
  itself instead.

- **`evalsolution` is removed.** Build a `Spline(basis, coeffs)` from the components or use the
  `Potential` functor.

- **A non-periodic basis without a boundary condition now fails loudly.** The branch that took one
  left two struct fields undefined and could not have worked. `PoissonSolverSpline` now checks the
  basis during construction and raises an `ArgumentError` naming the cause: a basis that represents
  the constants — the clamped basis, a Neumann recombination — has a singular stiffness matrix,
  because ``-\phi'' = \rho`` does not determine a constant. Only the periodic basis, whose constant
  mode the solver shifts out, and the Dirichlet basis, which has none, are accepted.

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

- **A periodic spline solver cannot be built at `Float32`.** `SimpleSplines.CirculantMass` checks
  that the assembled mass matrix is circulant against an absolute tolerance of `1e-10`, which
  `Float32` assembly noise exceeds, so construction fails with `ArgumentError: the mass matrix is
  not circulant to within 1.0e-10` before this package sees the matrix. The tolerance needs to
  scale with `eps(eltype)` upstream. The `Float32` Dirichlet spline and `Float32` grid solvers both
  work.
