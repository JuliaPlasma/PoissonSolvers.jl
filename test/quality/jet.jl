using JET
using PoissonSolvers
using Test

# Static optimisation analysis of the hot paths: every function of `src/` that a test file
# asserts with `@allocated`, at the concrete argument types those tests pass. A runtime
# dispatch on one of these paths is what the allocation tests measure only indirectly.
#
# Each further element type that a test outside `test/quality/` passes directly to the same
# method gets one line too, and so does each view at an offset it solves with: `Float32` for the
# `Potential` functor (test/potential.jl), and `Float32` and the unit-stride view for `solve!`
# (test/verification/tight_solves.jl). Another backend or another basis at the same element type
# needs no line. An element type that reaches a method only through another function gets none:
# the `Float32` solves of `solve` (test/matrixfree.jl) and of the `Potential` constructor
# (test/potential.jl) reach `solve!` that way.

if isdefined(JET, :JET_AVAILABLE) ? JET.JET_AVAILABLE :
   (isdefined(JET, :JET_LOADABLE) ? JET.JET_LOADABLE : true)
    m = (PoissonSolvers,)

    # test/poisson_fft.jl: solve! and the Potential functor on a grid basis
    fft = FFTWBasis((0.0, 1.0), 64)
    @test isempty(JET.get_reports(JET.report_opt(solve!,
        (Vector{Float64}, typeof(PoissonSolverFFT(fft)), Vector{Float64}); target_modules = m)))
    @test isempty(JET.get_reports(JET.report_opt(Potential(fft), (Float64,); target_modules = m)))

    # test/potential.jl: the Potential functor on a Float32 grid basis
    fft32 = FFTWBasis(Float32.((0.0, 1.0)), 64)
    @test isempty(JET.get_reports(JET.report_opt(Potential(fft32), (Float32,);
        target_modules = m)))

    # test/poisson_spline.jl: solve! on a periodic and a Dirichlet spline basis, which reach
    # different factorisations
    for b in (PeriodicBasisSpline((0.0, 1.0), 5, 32), DirichletBasisSpline((0.0, 1.0), 5, 32))
        @test isempty(JET.get_reports(JET.report_opt(solve!,
            (Vector{Float64}, typeof(PoissonSolverSpline(b)), Vector{Float64});
            target_modules = m)))
    end

    # test/matrixfree.jl: solve! of the matrix-free solver; the stencil order is a field, so
    # orders 2 and 4 give the same argument types
    grid = FiniteDifferenceBasis((0.0, 1.0), 64; order = 4)
    @test isempty(JET.get_reports(JET.report_opt(solve!,
        (Vector{Float64}, typeof(PoissonSolverMatrixFree(grid)), Vector{Float64});
        target_modules = m)))

    # test/verification/tight_solves.jl: solve! of each backend at `Float32`, and on the
    # unit-stride view at an offset the same file solves with, at both element types. Another
    # backend at the same element type is one of these methods; the periodic and the Dirichlet
    # spline are separate lines because they reach different factorisations.
    fft64 = PoissonSolverFFT(FFTWBasis((0.0, 1.0), 64))
    fft32 = PoissonSolverFFT(FFTWBasis(Float32.((0.0, 1.0)), 64))
    periodic64 = PoissonSolverSpline(PeriodicBasisSpline((0.0, 1.0), 5, 32))
    periodic32 = PoissonSolverSpline(PeriodicBasisSpline(Float32.((0.0, 1.0)), 5, 32))
    dirichlet64 = PoissonSolverSpline(DirichletBasisSpline((0.0, 1.0), 5, 32))
    dirichlet32 = PoissonSolverSpline(DirichletBasisSpline(Float32.((0.0, 1.0)), 5, 32))
    grid64 = PoissonSolverMatrixFree(grid)
    grid32 = PoissonSolverMatrixFree(FiniteDifferenceBasis(Float32.((0.0, 1.0)), 64; order = 4))

    view64 = typeof(view(zeros(66), 2:65))
    view32 = typeof(view(zeros(Float32, 66), 2:65))
    for p in (fft32, periodic32, dirichlet32, grid32)
        @test isempty(JET.get_reports(JET.report_opt(solve!,
            (Vector{Float32}, typeof(p), Vector{Float32}); target_modules = m)))
    end
    for p in (fft64, periodic64, dirichlet64, grid64)
        @test isempty(JET.get_reports(JET.report_opt(solve!,
            (Vector{Float64}, typeof(p), view64); target_modules = m)))
    end
    for p in (fft32, periodic32, dirichlet32, grid32)
        @test isempty(JET.get_reports(JET.report_opt(solve!,
            (Vector{Float32}, typeof(p), view32); target_modules = m)))
    end

    # The view as the *result* of a function right-hand side, which the same file's conversion
    # probe solves into as well. The function returns the other element type, as there.
    other(x::Float32) = 2.0 * x
    other(x::Float64) = Float32(2.0 * x)
    for p in (fft32, periodic32, dirichlet32, grid32)
        @test isempty(JET.get_reports(JET.report_opt(solve!,
            (view32, typeof(p), typeof(other)); target_modules = m)))
    end
    for p in (fft64, periodic64, dirichlet64, grid64)
        @test isempty(JET.get_reports(JET.report_opt(solve!,
            (view64, typeof(p), typeof(other)); target_modules = m)))
    end
else
    @test_skip "JET does not work on Julia $VERSION"  # aviatesk/JET.jl#681
end
