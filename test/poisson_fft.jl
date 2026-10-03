using PoissonSolvers
using PoissonSolvers: nearest_index, nearest_indices
using Random
using SimpleSplines
using Test

Random.seed!(1234)

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

    # CI runs the suite with --check-bounds=yes, which inflates allocations; the assertions are
    # therefore made only where bounds checking is at its default.
    if Base.JLOptions().check_bounds == 0
        @test solve_bytes == 0
        @test eval_bytes == 0
    end
end

@testset "a wrong length is a DimensionMismatch" begin
    # Every backend checks the three lengths before it reaches a plan or a factorisation, so a
    # wrong one is a DimensionMismatch naming the three counts.
    for T in (Float64, Float32)
        solver = PoissonSolverFFT(FFTWBasis(T.((0.0, 1.0)), 64))

        @test_throws DimensionMismatch solve!(zeros(T, 64), solver, zeros(T, 63))
        @test_throws DimensionMismatch solve!(zeros(T, 63), solver, zeros(T, 64))

        # The three lengths are distinct here, so the message pins which length is which.
        err = try
            solve!(zeros(T, 62), solver, zeros(T, 63))
        catch e
            e
        end
        @test err isa DimensionMismatch
        @test err.msg ==
              "the solver has 64 degrees of freedom, but the right-hand side has 63 and the result 62"

        # A zero-length result and an equal-but-wrong result disagree with the solver too. A check
        # that compared only the result with the right-hand side would accept both.
        @test_throws DimensionMismatch solve!(zeros(T, 0), solver, zeros(T, 0))
        @test_throws DimensionMismatch solve!(zeros(T, 63), solver, zeros(T, 63))

        # The plans are UNALIGNED so that they accept a strided view at an offset; a wrong-length
        # view must reach the length check rather than FFTW.
        long = zeros(T, 65)
        @test_throws DimensionMismatch solve!(zeros(T, 64), solver, view(long, 2:64))
        @test_throws DimensionMismatch solve!(view(long, 2:64), solver, zeros(T, 64))
    end
end

@testset "the function path reduces to the sampled vector path" begin
    # The function path reduces `f` through the shared hook, which for a grid samples the grid
    # points the vector path uses; the two are therefore `==`.
    for T in (Float64, Float32)
        basis = FFTWBasis(T.((0.0, 1.0)), 64)
        solver = PoissonSolverFFT(basis)
        f = x -> T(4π^2 * sin(2π * x) + 0.3 * cos(6π * x))
        rhs = f.(basis.xgrid[1:(end - 1)])

        @test solve(solver, f) == solve!(similar(rhs), solver, rhs)
        @test solve!(similar(rhs), solver, f) == solve!(similar(rhs), solver, rhs)

        # A strided view at an offset is a valid right-hand side too, because the plans accept one.
        wider = zeros(T, 65)
        wider[2:65] .= rhs
        @test solve!(similar(rhs), solver, view(wider, 2:65)) ==
              solve!(similar(rhs), solver, rhs)
    end
end

@testset "inference on both right-hand side kinds" begin
    function probe(::Type{T}) where {T}
        basis = FFTWBasis(T.((0.0, 1.0)), 64)
        solver = PoissonSolverFFT(basis)
        ρ = rand(T, length(solver))
        f = x -> T(4π^2 * sin(2π * x))
        (@inferred(solve!(similar(ρ), solver, ρ)),
            @inferred(solve!(similar(ρ), solver, f)),
            @inferred(solve(solver, ρ)), @inferred(solve(solver, f)))
    end
    for T in (Float64, Float32)
        a, b, c, d = probe(T)
        @test a isa Vector{T} && b isa Vector{T}
        @test c isa Vector{T} && d isa Vector{T}
    end
end
