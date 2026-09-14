using PoissonSolvers
using PoissonSolvers: nearest_index, nearest_indices
using SimpleSplines
using Test

@testset "spectral solve" begin
    # The transform diagonalises the periodic Laplacian exactly, so a single sine is reproduced
    # to round-off rather than to a discretisation error. Anything much larger means the symbol
    # or the normalisation is wrong.
    for (domain, f, u) in (((0.0, 1.0), x -> 4π^2 * sin(2π * x), x -> sin(2π * x)),
        ((-2π, 2π), sin, sin))
        basis = FFTWBasis(domain, 64)
        solver = PoissonSolverFFT(basis)
        xs = basis.xgrid[1:(end - 1)]

        @test PoissonSolver(basis) isa PoissonSolverFFT
        @test length(solver) == length(basis) == 64

        φ = solve(solver, f)
        @test φ ≈ u.(xs)
        @test solve(solver, f.(xs)) ≈ φ
        @test solve!(similar(φ), solver, f) ≈ φ

        # The constant mode is not determined by the equation, and this solver drops it.
        @test sum(φ) / length(φ) ≈ 0 atol=1e-12
    end
end

@testset "grid indexing" begin
    for domain in ((0.0, 1.0), (-2π, 2π))
        basis = FFTWBasis(domain, 64)

        @test nearest_index(basis, domain[begin]) == 1
        @test nearest_index(basis, domain[end]) == 1

        for x in range(domain[begin], domain[end], 17)[1:(end - 1)]
            i, j = nearest_indices(basis, x)
            @test basis.xgrid[i] ≤ x
            @test basis.xgrid[j] ≥ x
        end
    end
end

@testset "evaluation and derivative" begin
    basis = FFTWBasis((0.0, 1.0), 64)
    potential = Potential(basis, x -> 4π^2 * sin(2π * x))
    solution = potential.potential

    @test solution(0.25) == potential(0.25)
    @test solution(0.25, 0) == solution(0.25)

    # The derivative is the difference quotient of the two coefficients bracketing `x`. That is
    # the contract, and it is exact, so it is what gets pinned. No convergence rate is asserted
    # for it: the quotient approximates φ' at the midpoint of the cell holding `x` rather than at
    # `x` itself, so refining the grid moves the sample points around inside their cells and the
    # observed rate is not clean first order.
    c = coefficients(solution)
    for x in (0.23, 0.5, 0.77)
        i, j = nearest_indices(basis, x)
        @test solution(x, 1) == (c[j] - c[i]) / basis.Δx
    end

    # It is still a derivative: on a fine grid it tracks the analytic one to within the O(Δx)
    # such a scheme can offer.
    fine = Potential(FFTWBasis((0.0, 1.0), 256), x -> 4π^2 * sin(2π * x))
    @test maximum(abs, fine.(0.1:0.1:0.9, 1) .- 2π .* cos.(2π .* (0.1:0.1:0.9))) < 0.1

    @test derivative(solution).(0.1:0.1:0.9) == potential.(0.1:0.1:0.9, 1)
    @test_throws ArgumentError solution(0.25, 2)
end

@testset "type stability and allocations" begin
    function probe()
        basis = FFTWBasis((0.0, 1.0), 64)
        solver = PoissonSolverFFT(basis)
        ρ = rand(length(solver))
        φ = similar(ρ)
        solve!(φ, solver, ρ)
        potential = Potential(basis, ρ)
        (@inferred(solve!(φ, solver, ρ)), @allocated(solve!(φ, solver, ρ)),
            @inferred(potential(0.3)), @allocated(potential(0.3)))
    end
    result, solve_bytes, value, eval_bytes = probe()
    @test result isa Vector{Float64}
    @test value isa Float64

    # `Pkg.test()` forces --check-bounds=yes up to Julia 1.12, which inflates allocations; the
    # assertions are therefore made only where bounds checking is at its default.
    if Base.JLOptions().check_bounds == 0
        @test solve_bytes == 0
        @test eval_bytes == 0
    end
end
