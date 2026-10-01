using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "JET" include("quality/jet.jl")
    @safetestset "FFT Solvers" include("poisson_fft.jl")
    @safetestset "Spline Solvers" include("poisson_spline.jl")
    @safetestset "Matrix-Free Solvers" include("matrixfree.jl")
    @safetestset "Potential" include("potential.jl")
end
