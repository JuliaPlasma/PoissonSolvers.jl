using PoissonSolvers
using PoissonSolvers: nearest_index, nearest_indices
using Test

ngrid = 64


domain = (0.0, 1.0)

sol(x) = sin(2π * x)
rhs(x) = 4π^2 * sin(2π * x)

b = FFTWBasis(domain, ngrid)
s = PoissonSolverFFT(b)
ρ = rhs.(b.xgrid[1:end-1])
ϕ = sol.(b.xgrid[1:end-1])
φ = solve(s, ρ)
p = Potential(b, rhs)

@test s == PoissonSolver(b)

@test φ ≈ ϕ
@test φ == p.coefficients

@test nearest_index(b, domain[begin]) == 1
@test nearest_index(b, domain[end]) == 1

for x in (0.1, 0.25, 0.333333, 0.5)
    i,j = nearest_indices(b, x)
    @test b.xgrid[i] ≤ x
    @test b.xgrid[j] ≥ x
end


domain = (-2π, +2π)

sol(x) = sin(x)
rhs(x) = sin(x)

b = FFTWBasis(domain, ngrid)
s = PoissonSolverFFT(b)
ρ = rhs.(b.xgrid[1:end-1])
ϕ = sol.(b.xgrid[1:end-1])
φ = solve(s, ρ)
p = Potential(b, rhs)

@test s == PoissonSolver(b)

@test φ ≈ ϕ
@test φ == p.coefficients

@test nearest_index(b, domain[begin]) == 1
@test nearest_index(b, domain[end]) == 1

for x in (-4.0, -3.14, -0.5, 0.0, +0.5, +3.14, +4.0)
    i,j = nearest_indices(b, x)
    # println(" x = $x, i = $i, j = $j, xi = $(b.xgrid[i]), xj = $(b.xgrid[j])")
    @test b.xgrid[i] ≤ x
    @test b.xgrid[j] ≥ x
end
