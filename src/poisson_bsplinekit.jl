
using BSplineKit
using LinearAlgebra
using SparseArrays


struct PoissonSolverBSplineKit{DT<:Real,ST} <: PoissonSolver{DT}
    p::Int
    domain::Tuple{DT,DT}

    knots::Vector{DT}
    basis::ST

    M::LinearAlgebra.Hermitian{DT,SparseArrays.SparseMatrixCSC{DT,Int}}
    S::LinearAlgebra.Hermitian{DT,SparseArrays.SparseMatrixCSC{DT,Int}}

    P::Matrix{DT}
    R::Matrix{DT}

    Mfac::LinearAlgebra.Cholesky{DT,SparseArrays.SparseMatrixCSC{DT,Int}}
    Sfac::LinearAlgebra.Cholesky{DT,SparseArrays.Matrix{DT}}

    ρ::Vector{DT}
    ϕ::Vector{DT}
    rhs::Vector{DT}

    function PoissonSolverBSplineKit(domain::Tuple{DT,DT}, order::Int, n::Int) where {DT}
        knots = range(domain[begin], domain[end], n + 1)

        B = PeriodicBSplineBasis(BSplineOrder(order), knots)
        M = galerkin_matrix(B)
        S = galerkin_matrix(B, (Derivative(1), Derivative(1)))

        A = ones(n)
        R = A * A' / (A' * A)
        P = Matrix(I, n, n) .- R

        Mfac = cholesky!(M)
        Sfac = cholesky!(S .+ R)

        new{DT,typeof(B)}(order, domain, knots, B, M, S, P, R, Mfac, Sfac, zeros(DT, n), zeros(DT, n), zeros(DT, n))
    end
end


Base.length(p::PoissonSolverBSplineKit) = length(p.basis)


# compute RHS of L2 projection of the particle samples
function density!(p::PoissonSolverBSplineKit, points, weights)
    b = Splines.PeriodicVector(p.rhs)
    b .= 0

    for (x, w) in zip(points, weights)
        ilast, bs = p.basis(x)  # same as `evaluate_all`
        # Iterate over evaluated basis functions.
        # The indices of the evaluated basis functions are ilast:-1:(ilast - k + 1),
        # where k is the spline order.
        for (δi, bi) ∈ pairs(bs)
            i = ilast + 1 - δi
            b[i] += w * bi
        end
    end

    return p
end

function solve!(p::PoissonSolverBSplineKit, rhs::AbstractVector)
    ldiv!(p.ϕ, p.Sfac, rhs)
    return p
end

function solve!(p::PoissonSolverBSplineKit, x::AbstractVector, w::AbstractVector)
    density!(p, x, w)
    solve!(p, - p.P * p.rhs)
    ldiv!(p.ρ, p.Mfac, p.rhs)
end

function solve!(p::PoissonSolverBSplineKit, x::AbstractMatrix, w::AbstractMatrix)
    @assert size(x, 1) == 1
    @assert size(w, 1) == 1

    X = reshape(x, size(x, 2))
    W = reshape(w, size(w, 2))

    solve!(p, X, W)
end

function eval_density(p::PoissonSolverBSplineKit, x::Real)
    ρ = Spline(p.basis, p.ρ)
    ρ(x)
end

function eval_potential(p::PoissonSolverBSplineKit, x::Real)
    ϕ = Spline(p.basis, p.ϕ)
    ϕ(x)
end

function eval_field(p::PoissonSolverBSplineKit, x::Real)
    E = Derivative(1) * Spline(p.basis, p.ϕ)
    E(x)
end
