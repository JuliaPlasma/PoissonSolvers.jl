
using BSplineKit
using FFTW
using LinearAlgebra
using SparseArrays
using ToeplitzMatrices


function PeriodicBasisBSplineKit(domain, order, n)
    knots = range(domain[begin], domain[end], n + 1)
    PeriodicBSplineBasis(BSplineOrder(order), knots)
end


PoissonSolution(basis::AbstractBSplineBasis, coeffs::AbstractVector) = Spline(basis, coeffs)

function evalsolution(basis::AbstractBSplineBasis, coeffs::AbstractVector, x::Real)
    ϕ = Spline(basis, coeffs)
    ϕ(x)
end



struct PoissonSolverBSplineKit{DT<:Real, CT<:Complex, ST} <: PoissonSolver{DT}
    basis::ST

    M::Circulant{DT}
    S::Circulant{DT}

    Mfac::ToeplitzMatrices.ToeplitzFactorization{DT, Circulant{DT, SparseArrays.SparseVector{DT, Int}}, CT, FFTW.cFFTWPlan{CT, -1, true, 1, Tuple{Int}}}
    Sfac::ToeplitzMatrices.ToeplitzFactorization{DT, Circulant{DT, SparseArrays.SparseVector{DT, Int}}, CT, FFTW.cFFTWPlan{CT, -1, true, 1, Tuple{Int}}}

    P::Circulant{DT}
    R::Circulant{DT}

    function PoissonSolverBSplineKit(basis)
        M = galerkin_matrix(basis)
        S = galerkin_matrix(basis, (Derivative(1), Derivative(1)))

        Mcirc = Circulant(M[1,:])
        Scirc = Circulant(S[1,:])

        Mfac = factorize(Mcirc)

        println(typeof(Mfac))

        if typeof(basis) <: PeriodicBSplineBasis
            n = length(basis)
            A = ones(n)
            R = A * A' / (A' * A)
            P = Matrix(I, n, n) .- R
            
            Rcirc = Circulant(R[1,:])
            Pcirc = Circulant(P[1,:])

            Sfac = factorize(Scirc + Rcirc)

            new{eltype(M), complex(eltype(M)), typeof(basis)}(basis, Mcirc, Scirc, Mfac, Sfac, Pcirc, Rcirc)
        else
            Sfac = factorize(Scirc)

            new{eltype(M), complex(eltype(M)), typeof(basis)}(basis, Mcirc, Scirc, Mfac, Sfac)
        end
    end
end

PoissonSolver(basis::AbstractBSplineBasis) = PoissonSolverBSplineKit(basis)

knots(p::PoissonSolverBSplineKit) = BSplineKit.knots(p.basis)
order(p::PoissonSolverBSplineKit) = BSplineKit.order(p.basis)

Base.length(p::PoissonSolverBSplineKit) = length(p.basis)

isperiodic(b::AbstractBSplineBasis) = false
isperiodic(b::PeriodicBSplineBasis) = true
isperiodic(p::PoissonSolverBSplineKit) = isperiodic(p.basis)


function solve!(result::AbstractVector, p::PoissonSolverBSplineKit, rhs::AbstractVector)
    if isperiodic(p)
        # ldiv!(result, p.Sfac, p.P * rhs)
        result .= p.Sfac \ (p.P * rhs)
    else
        # ldiv!(result, p.Sfac, rhs)
        result .= p.Sfac \ rhs
    end
    return result
end

function solve!(result::AbstractVector, p::PoissonSolverBSplineKit, rhs::Base.Callable)
    solve!(result, p, galerkin_projection(rhs, p.basis))
end


function solve(p::PoissonSolverBSplineKit, rhs::AbstractVector)
    if isperiodic(p)
        return p.Sfac \ (p.P * rhs)
    else
        return p.Sfac \ rhs
    end
end

function solve(p::PoissonSolverBSplineKit, rhs::Base.Callable)
    solve(p, galerkin_projection(rhs, p.basis))
end
