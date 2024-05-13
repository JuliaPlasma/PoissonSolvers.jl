using SafeTestsets

@safetestset "FFT Solvers                                                                     " begin include("fft_tests.jl") end
@safetestset "BSplineKit Solvers                                                              " begin include("bsplinekit_tests.jl") end
