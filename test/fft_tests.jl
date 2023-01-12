using PoissonSolvers

ngrid = 64

domain = (-2π, +2π)

sol(x) = + sin(x)
rhs(x) = - sin(x)

p = PoissonSolverFFT(domain, ngrid)
ρ = rhs.(p.xgrid[1:end-1])
ϕ = sol.(p.xgrid[1:end-1])

solve!(p, ρ)

@test ϕ ≈ p.ϕ


domain = (0.0, 1.0)

sol(x) = + sin(2π * x)
rhs(x) = - 4π^2 * sin(2π * x)

p = PoissonSolverFFT(domain, ngrid)
ρ = rhs.(p.xgrid[1:end-1])
ϕ = sol.(p.xgrid[1:end-1])

solve!(p, ρ)

@test ϕ ≈ p.ϕ
