module PoissonSolvers

# Extended on this package's own solution and solver types, so that a caller who has `using
# PoissonSolvers` can reach a solution the same way they would reach a `SimpleSplines.Spline`.
import SimpleSplines: basis, coefficients, derivative

export basis, coefficients, derivative

include("poisson.jl")

export PoissonSolver
export solve!, solve

include("potential.jl")

export Potential

include("poisson_fft.jl")

export FFTWBasis
export PoissonSolverFFT

include("poisson_spline.jl")

export PeriodicBasisSpline, DirichletBasisSpline
export PoissonSolverSpline

include("matrixfree.jl")

export _apply_Δₓ!, _apply_Δₓ₄!, _apply_Rₓ!
export FiniteDifferenceBasis
export PoissonSolverMatrixFree

end
