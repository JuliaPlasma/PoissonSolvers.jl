
using BSplineKit
using LinearAlgebra
using SparseArrays


function PeriodicBasisBSplineKit(domain, order, n)
    knots = range(domain[begin], domain[end], n + 1)
    PeriodicBSplineBasis(BSplineOrder(order), knots)
end


PoissonSolution(basis::AbstractBSplineBasis, coeffs::AbstractVector) = Spline(basis, coeffs)

function evalsolution(basis::AbstractBSplineBasis, coeffs::AbstractVector, x::Real)
    ϕ = Spline(basis, coeffs)
    ϕ(x)
end



struct PoissonSolverBSplineKit{DT<:Real, ST} <: PoissonSolver{DT}
    basis::ST

    M::LinearAlgebra.Hermitian{DT,SparseArrays.SparseMatrixCSC{DT,Int}}
    S::LinearAlgebra.Hermitian{DT,SparseArrays.SparseMatrixCSC{DT,Int}}

    Mfac::LinearAlgebra.Cholesky{DT,SparseArrays.SparseMatrixCSC{DT,Int}}
    Sfac::LinearAlgebra.Cholesky{DT,SparseArrays.Matrix{DT}}

    P::Matrix{DT}
    R::Matrix{DT}

    function PoissonSolverBSplineKit(basis)
        M = galerkin_matrix(basis)
        S = galerkin_matrix(basis, (Derivative(1), Derivative(1)))

        Mfac = cholesky!(M)

        if typeof(basis) <: PeriodicBSplineBasis
            n = length(basis)
            A = ones(n)
            R = A * A' / (A' * A)
            P = Matrix(I, n, n) .- R

            Sfac = cholesky!(S .+ R)

            new{eltype(M), typeof(basis)}(basis, M, S, Mfac, Sfac, P, R)
        else
            Sfac = cholesky!(S)

            new{eltype(M), typeof(basis)}(basis, M, S, Mfac, Sfac)
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
        ldiv!(result, p.Sfac, p.P * rhs)
    else
        ldiv!(result, p.Sfac, rhs)
    end
    return result
end

function solve!(result::AbstractVector, p::PoissonSolverBSplineKit, rhs::Base.Callable)
    solve!(result, p, galerkin_projection(rhs, p.basis))
end


function solve(p::PoissonSolverBSplineKit, rhs::AbstractVector)
    if isperiodic(p)
        return p.Sfac \ p.P * rhs
    else
        return p.Sfac \ rhs
    end
end

function solve(p::PoissonSolverBSplineKit, rhs::Base.Callable)
    solve(p, galerkin_projection(rhs, p.basis))
end
