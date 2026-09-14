# PoissonSolvers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPlasma.github.io/PoissonSolvers.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaPlasma.github.io/PoissonSolvers.jl/latest/)
[![Build Status](https://github.com/JuliaPlasma/PoissonSolvers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaPlasma/PoissonSolvers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaPlasma/PoissonSolvers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaPlasma/PoissonSolvers.jl)

Solvers for the one-dimensional Poisson equation `-Δϕ = ρ`, with two backends:

- **a spectral solver** on a uniform periodic grid, which diagonalises the Laplacian by FFT and
  evaluates the solution at the nearest grid point;
- **a B-spline Galerkin solver** on the periodic and homogeneous-Dirichlet bases of
  [SimpleSplines](https://github.com/JuliaDEC/SimpleSplines.jl), of any order, which converges at
  the order of the basis and evaluates wherever it is asked.

A `Potential` pairs a solution with the solver that produced it, so a changing source is re-solved
in place. Every solve and every re-solve allocates nothing, on both backends.

```julia
using PoissonSolvers

ρ(x) = 4π^2 * sin(2π * x)

basis = PeriodicBasisSpline((0.0, 1.0), 5, 32)   # order k = 5, i.e. degree 4, on 32 cells
ϕ = Potential(basis, ρ)

ϕ(0.25)         # the potential
ϕ(0.25, 1)      # its first derivative

PoissonSolvers.update!(ϕ, ρ)   # re-solve in place for a new source
```

`update!` is deliberately not exported: the name is a common generic, and a caller that has
another one in scope should not have to disambiguate.

Use `DirichletBasisSpline` where the problem is not periodic, and `FFTWBasis` for the grid solver.


## Development

### Git hooks

Two hooks live in `.githooks`. They are **not active in a fresh clone** — `core.hooksPath` is local
configuration and does not travel with a push — so enable them once per clone:

```sh
git config core.hooksPath .githooks
```

**`pre-commit`** acts on **staged `.jl` files only**, and exits immediately when a commit stages
none, so a documentation- or workflow-only commit is not slowed down by it:

- **JuliaFormatter `--check`**, honouring this repository's own `.JuliaFormatter.toml` — **blocks**
  the commit. Formatting is mechanical and always fixable.
- **`fatou lint`**, when `fatou` is installed — **advisory only**, and deliberately so: its
  `unused-import` rule does not follow `include`, so it flags the load-bearing imports of every
  module file.
- **`using <Package>`**, which catches a syntax error or a broken `include` — **blocks**.

**`pre-push`** runs the full test suite with `--check-bounds=auto`, but **only when pushing to
`main` or `master`**; a topic branch is left to CI. It prints nothing for **10–30 minutes**, which
looks exactly like a network hang and is not one. If you do interrupt it, check for an orphaned
Julia process that the killed hook left behind.

Either hook can be bypassed for a single command with `--no-verify`, for a change you know it does
not apply to:

```sh
git commit --no-verify
git push --no-verify
```

The hooks are generated from one shared copy and are byte-identical across the related
repositories, so edit them there rather than here — a local edit is silently undone by the next
install.
