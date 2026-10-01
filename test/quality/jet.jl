using JET
using PoissonSolvers
using Test

# Static optimisation analysis of the hot paths: every function of `src/` that a test file
# asserts with `@allocated`, at the concrete argument types those tests pass. A runtime
# dispatch on one of these paths is what the allocation tests measure only indirectly.
#
# Each further element type that a test outside `test/quality/` passes directly to the same
# method gets one line too: `Float32` for the `Potential` functor (test/potential.jl). An element
# type that reaches a method only through another function gets none: the `Float32` solves of
# `solve` (test/matrixfree.jl) and of the `Potential` constructor (test/potential.jl) reach
# `solve!` that way.

if isdefined(JET, :JET_AVAILABLE) ? JET.JET_AVAILABLE : JET.JET_LOADABLE
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
else
    @test_skip "JET does not work on Julia $VERSION"  # aviatesk/JET.jl#681
end
