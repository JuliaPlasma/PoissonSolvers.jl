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

### K2 · A periodic spline solver cannot be built at `Float32`.

- **location:** —
- **evidence:** `SimpleSplines.CirculantMass` checks
  that the assembled mass matrix is circulant against an absolute tolerance of `1e-10`, which
  `Float32` assembly noise exceeds, so construction fails with `ArgumentError: the mass matrix is
  not circulant to within 1.0e-10` before this package sees the matrix. The tolerance needs to
  scale with `eps(eltype)` upstream. The `Float32` Dirichlet spline and `Float32` grid solvers both
  work.
- **kind:** upstream
- **found:** 2026-09-14

### K3 · Revise prints EMFILE errors in the test log

- **location:** `test/quality/jet.jl`
- **evidence:** JET 0.12 loads Revise, and its file watcher runs out of file handles. On Julia
  1.13.1, `grep -c 'UNHANDLED TASK ERROR'` counts 5 blocks in the `run-tests.jl` log of a full
  `Pkg.test()` of the branch that adds `test/quality/jet.jl`, each an
  `IOError: FolderMonitor: too many open files (EMFILE)` stack trace. The test environment of
  `origin/main` holds no JET, so no Revise. These blocks are not test failures, and the test totals
  do not change.
- **kind:** upstream
- **found:** 2026-10-01
