using SafeTestsets

@safetestset "Aqua                                                                            " begin
    include("aqua_tests.jl")
end
@safetestset "FFT Solvers                                                                     " begin
    include("fft_tests.jl")
end
@safetestset "Spline Solvers                                                                  " begin
    include("spline_tests.jl")
end
@safetestset "Matrix-Free Solvers                                                             " begin
    include("matrixfree_tests.jl")
end
@safetestset "Potential                                                                       " begin
    include("potential_tests.jl")
end
