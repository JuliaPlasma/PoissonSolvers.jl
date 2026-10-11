# Known issues

### K1 · Evaluating a spline potential at a scalar point allocates 96 bytes.

- **location:** `scripts/measure_allocations.jl`
- **evidence:** This is not in this
  package: `SimpleSplines.evaluate(basis, coefficients, x)` takes one `local_width` buffer per
  scalar call, because the degree is a struct field rather than a type parameter. The fix belongs
  upstream: `SimpleSplines.evaluate_all!`, which takes a caller-supplied buffer, allocates nothing.
  See `scripts/measure_allocations.jl` for both figures side by side. The grid backend evaluates
  with no allocation.
- **kind:** upstream
- **found:** 2026-09-14

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
