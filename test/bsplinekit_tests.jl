using BSplineKit
using PoissonSolvers

nknot = 32
order = 5


## test with manufactured solutions

domain = (-2π, +2π)
x = domain[begin]:0.1:domain[end]

sol(x) = sin(x)
rhs(x) = sin(x)

p = PoissonSolverBSplineKit(domain, order, nknot)
ρ = galerkin_projection(rhs, p.basis)
ϕ = sol.(x)

solve!(p, ρ)

@test ϕ ≈ eval_potential(p, x)  atol = 1e-3


domain = (0.0, 1.0)
x = domain[begin]:0.1:domain[end]

sol(x) = sin(2π * x)
rhs(x) = 4π^2 * sin(2π * x)

p = PoissonSolverBSplineKit(domain, order, nknot)
ρ = galerkin_projection(rhs, p.basis)
ϕ = sol.(x)

solve!(p, ρ)

@test ϕ ≈ eval_potential(p, x)  atol = 1e-3


## compare with legacy solver

nknot = 16
order = 3

domain = (0.0, 1.0)

ps = PoissonSolverBSplineKit(domain, order, nknot)
pr = PoissonSolverPBSplines(domain, order, nknot)

@test length(ps) == length(pr)


using AdaptiveRejectionSampling

n = 10000
μ = 0.0 
σ = 2.0
f = x -> exp(-0.5 * (4π * x - μ - 2π)^2 / σ^2) / sqrt(2π * σ^2)
x = run_sampler!(RejectionSampler(f, domain, max_segments=5), n)
w = ones(n) ./ n

solve!(ps, x, w)
solve!(pr, x, w)

x = domain[begin]:0.1:domain[end]

ϕ1 = eval_potential(ps, x)
ϕ2 = eval_potential(pr, x)

@test ϕ1 ≈ ϕ2  atol = 1E-4

E1 = eval_field(ps, x)
E2 = eval_field(pr, x)

@test E1 ≈ E2  atol = 1E-2
