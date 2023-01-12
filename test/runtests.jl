using PoissonSolvers
using Test

@testset "PoissonSolvers.jl" begin
    @testset "FFT Solvers" begin
        include("fft_tests.jl")
    end
    @testset "BSplineKit Solvers" begin
        include("bsplinekit_tests.jl")
    end
end
