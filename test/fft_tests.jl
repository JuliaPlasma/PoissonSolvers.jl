using PoissonSolvers

ngrid = 64

domain = (-2π, +2π)


sol(x) = + sin(x)
rhs(x) = - sin(x)

b = FFTWBasis(domain, ngrid)
s = PoissonSolverFFT(b)
ρ = rhs.(b.xgrid[1:end-1])
ϕ = sol.(b.xgrid[1:end-1])
φ = solve(s, ρ)
p = Potential(b, rhs)

@test s == PoissonSolver(b)

@test φ ≈ ϕ
@test φ == p.coefficients


domain = (0.0, 1.0)

sol(x) = + sin(2π * x)
rhs(x) = - 4π^2 * sin(2π * x)

b = FFTWBasis(domain, ngrid)
s = PoissonSolverFFT(b)
ρ = rhs.(b.xgrid[1:end-1])
ϕ = sol.(b.xgrid[1:end-1])
φ = solve(s, ρ)
p = Potential(b, rhs)

@test s == PoissonSolver(b)

@test φ ≈ ϕ
@test φ == p.coefficients
