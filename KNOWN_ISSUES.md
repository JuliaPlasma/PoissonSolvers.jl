# Known issues

### K1 · Evaluating a spline potential at a scalar point allocates one buffer per call.

- **location:** `src/potential.jl`
- **evidence:** `p(x)`, `p(x, d)` and `derivative(p, d)(x)` on a spline basis call SimpleSplines'
  `Spline`, whose `evaluate` takes one `local_width` buffer per scalar call: 96 B at `Float64` and
  80 B at `Float32` at order 5 with SimpleSplines 0.4.0; the size grows with the order.
  `scripts/measure_allocations.jl` gives the `Float64` figure beside the grid backend's, which
  evaluates with no allocation. SimpleSplines 0.4's `evaluate!` takes a caller's buffer and
  allocates nothing, so a buffer held by `Potential` would remove the allocation. The author
  accepts it: point evaluation of a `Potential` serves mainly plotting, and the hot paths of the
  callers evaluate through `evaluate_all!` with their own buffer.
- **kind:** defect
- **found:** 2026-09-14; accepted by the author on 2026-11-03

### K2 · A `Float32` periodic spline is rejected at a fine resolution.

- **location:** —
- **evidence:** `SimpleSplines`' `_reciprocal_eigenvalues` treats an eigenvalue below
  `n * eps(T) * maximum(abs, ĉ)` as zero, which is `pinv`'s own relative cutoff, and the smallest
  nonzero eigenvalue of a periodic stiffness matrix falls under it at a fine resolution. A
  `PeriodicBasisSpline` at `Float32` builds at n = 256 for orders 2 to 5, and fails with
  `ArgumentError: the matrix is singular beyond the constants: a second mode of the transform
  vanishes as well` at n = 512 at order 2 and at n = 1024 at order 4. The `Float32` Dirichlet
  spline and the `Float32` grid solvers build at n = 256, 512 and 1024. The fix is upstream.
- **kind:** upstream
- **found:** 2026-09-14

### K3 · Revise prints EMFILE errors in the test log.

- **location:** `test/quality/jet.jl`
- **evidence:** JET 0.12 loads Revise, and its file watcher runs out of file handles. On Julia
  1.13.1, `grep -c 'UNHANDLED TASK ERROR'` counts 5 blocks in the `run-tests.jl` log of a full
  `Pkg.test()` of the branch that adds `test/quality/jet.jl`, each an
  `IOError: FolderMonitor: too many open files (EMFILE)` stack trace. The same count on a full run
  of `origin/main` at eda0615 gives 0: its test environment holds no JET, so no Revise. These
  blocks are not test failures, and the test totals do not change.
- **kind:** upstream
- **found:** 2026-10-01

### K4 · A constant right-hand side near `floatmax` solves to `NaN`, not to zero.

- **location:** `src/poisson_fft.jl:139`
- **evidence:** `solve(PoissonSolver(FFTWBasis(T.((0, 1)), 64)), x -> floatmax(T) / 2)` returns a
  vector that is not all finite, in `Float32` and in `Float64`. A constant is in the kernel of the
  periodic operator and the k = 0 factor is zero instead of `1/k²`, so the mean-free solution is
  zero; the transform's range instead overflows the DC mode and `Inf * 0` reaches the result. The
  input is outside the density and scale domain of `test/verification/tight_solves.jl` (`1`, `3e-8`,
  `3e8`).
- **kind:** found late
- **found:** 2026-10-11; reproduced on base commit `6146203` as well as this branch.

### K5 · A periodic spline solver does not build beside FFTA 0.3.2.

- **location:** —
- **evidence:** `PeriodicBasisSpline` builds its mass operator through SimpleSplines'
  `CirculantMass`, which calls `plan_rfft(c; flags)`. FFTA 0.3.2 adds a `plan_rfft` method that is
  ambiguous with FFTW 1.10's, so with FFTA loaded, as CairoMakie loads it, the build raises
  `MethodError: … plan_rfft … is ambiguous` with SimpleSplines 0.3.4. The FFT solver has the same
  call and builds beside FFTA; `test/integration/ffta.jl` checks that and runs last in `core`, so
  that the spline tests do not run beside FFTA. The fix is upstream, in JuliaDEC/SimpleSplines.jl
  PR #43.
- **kind:** upstream
- **found:** 2026-11-02
