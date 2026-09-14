using PoissonSolvers
using SimpleSplines
using Test

# The spline solvers are verified by convergence rather than by a single tolerance. A B-spline
# Galerkin discretisation of order k converges at order k, so the observed rate is what says the
# assembly, the boundary treatment and the mass solve are all right; one small number at one
# resolution would pass just as happily with a subtly wrong load vector.
const XT = range(0, 1, 201)

"""Maximum error of the solver on `basis` against the exact solution `u`, right-hand side `f`."""
function solution_error(basis, f, u)
    s = Spline(basis, solve(PoissonSolverSpline(basis), f))
    maximum(abs, s.(XT) .- u.(XT))
end

function convergence_rates(mkbasis, order, f, u)
    errs = [solution_error(mkbasis((0.0, 1.0), order, n), f, u) for n in (8, 16, 32, 64)]
    [log2(errs[i] / errs[i + 1]) for i in 1:(length(errs) - 1)]
end

@testset "periodic solver" begin
    f(x) = 4π^2 * sin(2π * x)
    u(x) = sin(2π * x)

    for order in (3, 5)
        @test all(≥(order - 0.5), convergence_rates(PeriodicBasisSpline, order, f, u))
    end

    basis = PeriodicBasisSpline((0.0, 1.0), 5, 32)
    solver = PoissonSolverSpline(basis)

    @test PoissonSolver(basis) isa PoissonSolverSpline
    @test length(solver) == nbasis(basis) == 32
    @test PoissonSolvers.basis(solver) === basis

    # The solution of a periodic Poisson problem is fixed only up to a constant, and this solver
    # picks the mean-free one. Nothing else in the suite pins that choice.
    φ = solve(solver, f)
    @test sum(φ) / length(φ) ≈ 0 atol=1e-12

    @test solve!(similar(φ), solver, f) ≈ φ
    @test solve!(similar(φ), solver, PoissonSolvers.loadvector(solver, f)) ≈ φ

    # A right-hand side with a constant part has that part discarded, not diffused into the
    # answer: adding one changes nothing.
    ρ = PoissonSolvers.loadvector(solver, f)
    @test solve(solver, ρ .+ 1.0) ≈ φ
end

@testset "a second domain" begin
    # -φ" = sin on (-2π, 2π), which is the same equation with a different period, and catches a
    # solver that silently assumes the unit interval.
    basis = PeriodicBasisSpline((-2π, 2π), 5, 32)
    s = Spline(basis, solve(PoissonSolverSpline(basis), sin))
    xs = range(-2π, 2π, 201)
    @test maximum(abs, s.(xs) .- sin.(xs)) < 1e-3
end

@testset "Dirichlet solver" begin
    f(x) = π^2 * sin(π * x)
    u(x) = sin(π * x)

    for order in (3, 5)
        @test all(≥(order - 0.5), convergence_rates(DirichletBasisSpline, order, f, u))
    end

    basis = DirichletBasisSpline((0.0, 1.0), 5, 32)
    solver = PoissonSolverSpline(basis)
    s = Spline(basis, solve(solver, f))

    # Recombination is what enforces the boundary condition, so it holds exactly rather than
    # approximately — every basis function vanishes at both ends.
    @test s(0.0) == 0.0
    @test s(1.0) == 0.0

    # Two degrees of freedom fewer than the clamped basis of the same order and mesh.
    @test nbasis(basis) == nbasis(BSplineBasis(UniformMesh(32, (0.0, 1.0)), 4)) - 2

    # The stiffness matrix is non-singular here, so no mean is removed: a constant added to the
    # right-hand side does change the answer, unlike the periodic case above.
    ρ = PoissonSolvers.loadvector(solver, f)
    @test !isapprox(solve(solver, ρ .+ 1.0), solve(solver, ρ))
end

@testset "order and degree" begin
    # `order` is k = p + 1 at this package's boundary, while SimpleSplines takes the degree p.
    # An off-by-one here would still converge, just at the wrong rate, so pin it directly.
    for order in (2, 3, 5)
        @test degree(PeriodicBasisSpline((0.0, 1.0), order, 16)) == order - 1
        @test SimpleSplines.order(PeriodicBasisSpline((0.0, 1.0), order, 16)) == order
    end
end

@testset "type stability and allocations" begin
    function probe()
        basis = PeriodicBasisSpline((0.0, 1.0), 5, 32)
        solver = PoissonSolverSpline(basis)
        ρ = rand(length(solver))
        φ = similar(ρ)
        solve!(φ, solver, ρ)
        (@inferred(solve!(φ, solver, ρ)), @allocated(solve!(φ, solver, ρ)))
    end
    result, allocated = probe()
    @test result isa Vector{Float64}

    # `Pkg.test()` forces --check-bounds=yes up to Julia 1.12, which inflates allocations; the
    # assertion is therefore made only where bounds checking is at its default.
    if Base.JLOptions().check_bounds == 0
        @test allocated == 0
    end
end

@testset "rejected input" begin
    solver = PoissonSolverSpline(PeriodicBasisSpline((0.0, 1.0), 5, 32))
    @test_throws DimensionMismatch solve!(zeros(8), solver, rand(32))

    # A periodic basis of degree p needs more than p cells for the wrap to be well defined.
    @test_throws ArgumentError PeriodicBasisSpline((0.0, 1.0), 5, 3)

    # A basis that represents the constants has a singular stiffness matrix, and the solver says
    # so itself rather than letting the factorisation report it against the mass matrix.
    mesh = UniformMesh(32, (0.0, 1.0))
    @test_throws ArgumentError PoissonSolverSpline(BSplineBasis(mesh, 4))
    @test_throws ArgumentError PoissonSolverSpline(BSplineBasis(mesh, 4, Neumann()))
end

@testset "the stiffness matrix reaches the operator unmodified" begin
    # The singularity is deflated in the solve, so the assembly is what the operator is handed
    # and what it keeps. The alternative cure — a shift by the rank-one mean projector 𝟙𝟙ᵀ/n —
    # is a scalar added to every entry, so it makes a banded assembly structurally full and it
    # promotes a narrower element type to the scalar's own.
    b = PeriodicBasisSpline((0.0, 1.0), 4, 32)
    solver = PoissonSolverSpline(b)
    S = stiffness_matrix(SplineQuadrature(b))

    stored = count(!iszero, Matrix(SimpleSplines.mass_matrix(solver.stiffness)))
    @test stored == count(!iszero, Matrix(S))
    @test stored < nbasis(b)^2 ÷ 4
    @test eltype(solver.stiffness) == Float64

    # and the deflation is what makes that possible: the solution is still the mean-free one
    @test sum(solve(solver, randn(nbasis(b)))) ≈ 0 atol = 1e-12
end
