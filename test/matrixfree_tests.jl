using PoissonSolvers
using PoissonSolvers: _apply_L!, _apply_Δₓ!, _apply_Δₓ₄!, _apply_Rₓ!
using LinearAlgebra
using Test

source(x) = 4π^2 * sin(2π * x)
exact(x) = sin(2π * x)

# A smooth source with every Fourier mode present, so that no test passes because the grid
# happens to resolve a single mode exactly.
smooth(x) = exp(sin(2π * x)) + cos(4π * x)

gridpoints(b) = b.xgrid[1:(end - 1)]

# The observed orders of convergence of `errors`, each on a grid twice as fine as the last.
rates(errors) = [log2(errors[i] / errors[i + 1]) for i in 1:(length(errors) - 1)]

@testset "construction" begin
    for order in (2, 4)
        b = FiniteDifferenceBasis((0.0, 1.0), 64; order)
        solver = PoissonSolver(b)

        @test solver isa PoissonSolverMatrixFree
        @test basis(solver) === b
        @test length(solver) == length(b) == 64
        @test b.order == order
    end
    @test FiniteDifferenceBasis((0.0, 1.0), 64).order == 2
    @test_throws ArgumentError FiniteDifferenceBasis((0.0, 1.0), 64; order = 3)
    g = FFTWBasis((0.0, 1.0), 64)
    @test_throws ArgumentError FiniteDifferenceBasis(g.domain, g.xgrid, g.Δx, 3)
end

@testset "stencil order" begin
    # The stencils on their own, before any solve: each approximates φ'' to its nominal order.
    for (apply!, order) in ((_apply_Δₓ!, 2), (_apply_Δₓ₄!, 4))
        errors = map((16, 32, 64, 128)) do n
            b = FiniteDifferenceBasis((0.0, 1.0), n)
            y = similar(gridpoints(b))
            apply!(y, exact.(gridpoints(b)), b.Δx)
            maximum(abs, y .+ source.(gridpoints(b)))
        end
        @test all(r -> abs(r - order) < 0.1, rates(errors))
    end

    # `_apply_Lₓ₄!` fuses the fourth-order stencil with the projection, so it applies the same
    # operator as the solver composes from the separate stencils.
    solver = PoissonSolver(FiniteDifferenceBasis((0.0, 1.0), 64; order = 4))
    x = randn(64)
    y = similar(x)
    PoissonSolvers._apply_Lₓ₄!(y, x, solver.basis.Δx)
    @test y ≈ _apply_L!(similar(x), x, solver)
end

@testset "refinement sweep" begin
    for order in (2, 4)
        errors = map((16, 32, 64, 128)) do n
            b = FiniteDifferenceBasis((0.0, 1.0), n; order)
            maximum(abs, solve(PoissonSolver(b), source) .- exact.(gridpoints(b)))
        end
        @test all(r -> abs(r - order) < 0.1, rates(errors))
    end
end

@testset "the conjugate gradients solve the discrete system" begin
    # Independent of the discretisation error: the solver inverts its own operator, on a
    # mean-free vector with no smoothness at all, to round-off.
    for order in (2, 4), n in (64, 1024)

        solver = PoissonSolver(FiniteDifferenceBasis((0.0, 1.0), n; order))
        φ = randn(n)
        φ .-= sum(φ) / n
        ρ = _apply_L!(similar(φ), φ, solver)
        @test solve(solver, ρ) ≈ φ rtol=1e-10
    end
end

@testset "Float32 accuracy" begin
    # `reltol` bounds the residual, not the error, so this checks the error against the analytic
    # solution. exp(sin(2πx)) has every Fourier mode, so conjugate gradients cannot finish in a few
    # steps on a few eigenvectors, and the solver returns its mean-free part.
    potential(x) = exp(sin(2π * x))
    density(x) = 4π^2 * exp(sin(2π * x)) * (sin(2π * x) - cos(2π * x)^2)
    b = FiniteDifferenceBasis(Float32.((0.0, 1.0)), 256; order = 4)
    x = Float64.(gridpoints(b))
    reference = potential.(x) .- sum(potential.(x)) / length(x)
    φ = solve(PoissonSolver(b), Float32.(density.(x)))

    @test eltype(φ) == Float32
    @test norm(φ .- reference) / norm(reference) ≤ 1e-3
end

@testset "the same periodic problem as the FFT backend" begin
    # `smooth` has period 1, so each domain spans a whole number of periods.
    for domain in ((0.0, 1.0), (-1.0, 1.0))
        fd = PoissonSolver(FiniteDifferenceBasis(domain, 256; order = 4))
        fft = PoissonSolver(FFTWBasis(domain, 256))

        @test maximum(abs, solve(fd, smooth) .- solve(fft, smooth)) ≤ sqrt(eps())
    end
end

@testset "the nullspace projection matches the kernel deflation" begin
    b = FiniteDifferenceBasis((0.0, 1.0), 64; order = 4)
    fft = PoissonSolverFFT(FFTWBasis((0.0, 1.0), 64))
    spline = PoissonSolver(PeriodicBasisSpline((0.0, 1.0), 5, 32))
    ρ = smooth.(gridpoints(b))

    # `1 - R` is exactly what the FFT backend does to the right-hand side, where the factor of
    # the k = 0 mode is zero.
    Rρ = similar(ρ)
    _apply_Rₓ!(Rρ, ρ)
    @test all(==(sum(ρ) / length(ρ)), Rρ)
    deflated = similar(ρ)
    mul!(fft.ρ̂, fft.plan, ρ)
    fft.ρ̂[1] = 0
    mul!(deflated, fft.iplan, fft.ρ̂)
    @test ρ .- Rρ ≈ deflated rtol=1e-14

    # Every backend drops a constant in the source rather than failing on it, and returns the
    # mean-free solution.
    for solver in (PoissonSolver(b), fft, spline)
        @test solve(solver, x -> 1 + smooth(x)) ≈ solve(solver, smooth) rtol=1e-12
        @test maximum(abs, solve(solver, x -> 1.0)) ≤ 1e-12
    end
    @test sum(solve(PoissonSolver(b), smooth)) / length(b) ≈ 0 atol=1e-14
end

@testset "a wrong length is a DimensionMismatch" begin
    solver = PoissonSolver(FiniteDifferenceBasis((0.0, 1.0), 64))
    @test_throws DimensionMismatch solve!(zeros(64), solver, zeros(63))
    @test_throws DimensionMismatch solve!(zeros(63), solver, zeros(64))
end

@testset "non-convergence throws" begin
    b = FiniteDifferenceBasis((0.0, 1.0), 64)
    @test_throws ErrorException solve(PoissonSolverMatrixFree(b; maxiter = 2), smooth)
end

@testset "type stability and allocations" begin
    function probe(order)
        b = FiniteDifferenceBasis((0.0, 1.0), 64; order)
        solver = PoissonSolverMatrixFree(b)
        ρ = rand(length(solver))
        φ = similar(ρ)
        solve!(φ, solver, ρ)
        (@inferred(solve!(φ, solver, ρ)), @allocated(solve!(φ, solver, ρ)))
    end
    for order in (2, 4)
        result, bytes = probe(order)
        @test result isa Vector{Float64}

        # CI runs the suite with --check-bounds=yes, which inflates allocations; the assertion is
        # therefore made only where bounds checking is at its default.
        if Base.JLOptions().check_bounds == 0
            @test bytes == 0
        end
    end
end
