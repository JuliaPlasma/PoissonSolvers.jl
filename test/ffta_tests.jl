using PoissonSolvers
using FFTA: FFTA
using Test

# FFTA, loaded here as Makie loads it, adds `plan_rfft` and `plan_brfft` methods that make the
# plain FFTW calls ambiguous. This file runs last, so that no other file runs beside FFTA.
@testset "the FFT solver builds and solves beside FFTA" begin
    for DT in (Float32, Float64), n in (63, 64)

        basis = FFTWBasis((DT(0), DT(1)), n)
        solver = PoissonSolverFFT(basis)
        @test solver isa PoissonSolverFFT{DT}
        xs = basis.xgrid[1:(end - 1)]
        φ = solve(solver, DT(4π^2) .* sinpi.(2 .* xs))
        @test eltype(φ) === DT
        @test φ ≈ sinpi.(2 .* xs) rtol = 64 * eps(DT)
    end
end
