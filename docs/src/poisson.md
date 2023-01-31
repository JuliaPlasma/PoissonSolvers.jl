```@meta
CurrentModule = PoissonSolvers
```

# Poisson Solvers

Solve
```math
- \Delta \phi(x) = f(x)
```


```@example 1
using Plots
using PoissonSolvers

sol(x) = sin(π*x) * (1-x^2)
der(x) = π * cos(π*x) * (1-x^2) - 2 * x * sin(π*x)
rhs(x) = 4π * x * cos(π*x) + 2 * sin(π*x) + π^2 * (1-x^2) * sin(π*x)

domain = (-1., +1.)

x = LinRange(domain[begin], domain[end], 100)
```


## FFT Solver

```@example 1
b = FFTWBasis(domain, 256)
p = Potential(b, rhs)
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ(x)")
plot!(x, p.(x); xlims = domain, label = "Solution")
plot!(x, sol.(x); xlims = domain, label = "Reference")
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ'(x)")
plot!(x, p.(x, Derivative(1)); xlims = domain, label = "Derivative")
plot!(x, der.(x); xlims = domain, label = "Reference")
```


## B-Spline Solver

```@example 1
b = PeriodicBasisBSplineKit(domain, 3, 32)
p = Potential(b, rhs)
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ(x)")
plot!(x, p.(x); xlims = domain, label = "Solution")
plot!(x, sol.(x); xlims = domain, label = "Reference")
```

```@example 1
plot(xlabel = "x", ylabel = "ϕ'(x)")
plot!(x, p.(x, Derivative(1)); xlims = domain, label = "Derivative")
plot!(x, der.(x); xlims = domain, label = "Reference")
```
