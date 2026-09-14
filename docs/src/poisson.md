```@meta
CurrentModule = PoissonSolvers
```

# Poisson Solvers

Solve
```math
- \Delta \phi(x) = f(x)
```
in one dimension, with either a spectral solver on a uniform periodic grid or a B-spline Galerkin
solver. A [`Potential`](@ref) pairs the solution with the solver that produced it, so that it can
be re-solved in place as the source changes.

Evaluate a potential with `p(x)`, and its `d`-th derivative with `p(x, d)`.

```@example 1
using Plots
using PoissonSolvers

sol(x) = sin(2π * x)
der(x) = 2π * cos(2π * x)
rhs(x) = 4π^2 * sin(2π * x)

domain = (0.0, 1.0)

x = LinRange(domain[begin], domain[end], 200)
```

## FFT Solver

The grid solver is spectral in the coefficients, but evaluates at the nearest grid point, so the
plotted curve is a staircase of width `Δx`.

```@example 1
b = FFTWBasis(domain, 64)
p = Potential(b, rhs)
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ(x)")
plot!(x, p.(x); xlims = domain, label = "Solution")
plot!(x, sol.(x); xlims = domain, label = "Reference")
```

## Periodic B-Spline Solver

A spline of order ``k`` converges at order ``k``, and evaluates wherever it is asked.

```@example 1
b = PeriodicBasisSpline(domain, 5, 32)
p = Potential(b, rhs)
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ(x)")
plot!(x, p.(x); xlims = domain, label = "Solution")
plot!(x, sol.(x); xlims = domain, label = "Reference")
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ'(x)")
plot!(x, p.(x, 1); xlims = domain, label = "Derivative")
plot!(x, der.(x); xlims = domain, label = "Reference")
```

## Dirichlet B-Spline Solver

A periodic basis needs a periodic problem. For a source that does not repeat, use a basis
recombined to vanish at both ends, which solves the same equation with ``\phi(a) = \phi(b) = 0``.

```@example 1
dsol(x) = sin(π * x)
drhs(x) = π^2 * sin(π * x)

b = DirichletBasisSpline(domain, 5, 32)
p = Potential(b, drhs)
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ(x)")
plot!(x, p.(x); xlims = domain, label = "Solution")
plot!(x, dsol.(x); xlims = domain, label = "Reference")
```

The boundary condition is built into the basis rather than imposed on the solution, so it holds
exactly:

```@example 1
p(domain[begin]), p(domain[end])
```
