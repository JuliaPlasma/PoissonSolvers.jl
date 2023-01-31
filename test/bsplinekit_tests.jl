using BSplineKit
using PoissonSolvers
using PoissonSolvers: evalsolution

nknot = 32
order = 5


## test with manufactured solutions

domain = (-2π, +2π)
x = domain[begin]:0.1:domain[end]

sol(x) = sin(x)
rhs(x) = sin(x)

b = PeriodicBasisBSplineKit(domain, order, nknot)
s = PoissonSolverBSplineKit(b)
ϕ = sol.(x)
φ = solve(s, rhs)
p = Potential(b, rhs)

# @test s == PoissonSolver(b)

@test evalsolution(b, φ, x) ≈ ϕ  atol = 1e-3
@test evalsolution(b, φ, x) == p(x)
@test φ == p.coefficients



domain = (0.0, 1.0)
x = domain[begin]:0.1:domain[end]

sol(x) = sin(2π * x)
rhs(x) = 4π^2 * sin(2π * x)

b = PeriodicBasisBSplineKit(domain, order, nknot)
s = PoissonSolverBSplineKit(b)
ϕ = sol.(x)
φ = solve(s, rhs)
p = Potential(b, rhs)

# @test s == PoissonSolver(b)

@test evalsolution(b, φ, x) ≈ ϕ  atol = 1e-3
@test evalsolution(b, φ, x) == p(x)
@test φ == p.coefficients
