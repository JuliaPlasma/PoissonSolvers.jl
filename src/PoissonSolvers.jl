module PoissonSolvers

import BSplineKit: Derivative


export Derivative


include("poisson.jl")

export PoissonSolver
export solve!, solve


include("potential.jl")

export Potential


include("poisson_fft.jl")

export FFTWBasis
export PoissonSolverFFT, PeriodicBasisFFT


include("poisson_bsplinekit.jl")

export PoissonSolverBSplineKit, PeriodicBasisBSplineKit

end
