```@meta
CurrentModule = PoissonSolvers
```

# Poisson Solvers


```@example 1
using Distributions
using Plots
using PoissonSolvers
using Random
using StatsBase

np = 10000
xp = zeros(np)
wp = ones(np) ./ np
rand!(MersenneTwister(0), Normal(0.5, 0.1), xp)

x = LinRange(0, 1, 100)
y = fit(Histogram, xp, x).weights ./ length(x)
x = x[1:end-1] .+ (x[2] - x[1]) / 2

plot(x, y; xlims = (0,1), xlabel = "x", ylabel = "n", legend = :none)
```


## FFT Solver

```@example 1
p = PoissonSolverFFT((0.0, 1.0), 32)
solve!(p, xp, wp);
```


```@example 1
plot(x, eval_density(p, x); xlims = (0,1), xlabel = "x", ylabel = "ρ(x)", legend = :none)
```

```@example 1
plot(x, eval_potential(p, x); xlims = (0,1), xlabel = "x", ylabel = "ϕ(x)", legend = :none)
```

```@example 1
plot(x, eval_field(p, x); xlims = (0,1), xlabel = "x", ylabel = "E(x)", legend = :none)
```


## B-Spline Solver

```@example 1
p = PoissonSolverPBSplines((0.0, 1.0), 3, 32)
solve!(p, xp, wp);
```


```@example 1
plot(x, eval_density(p, x); xlims = (0,1), xlabel = "x", ylabel = "ρ(x)", legend = :none)
```

```@example 1
plot(x, eval_potential(p, x); xlims = (0,1), xlabel = "x", ylabel = "ϕ(x)", legend = :none)
```

```@example 1
plot(x, eval_field(p, x); xlims = (0,1), xlabel = "x", ylabel = "E(x)", legend = :none)
```
