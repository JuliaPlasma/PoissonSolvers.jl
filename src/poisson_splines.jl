
using LinearAlgebra

function massmatrix(S::PBSpline{T}, q=S.p + 1) where {T}
    M = zeros(T, S.nₕ, S.nₕ)
    for i in 1:S.nₕ
        for k in 0:S.p
            M[1, i] += integrate_gausslegendre(x -> (eval_PBSpline.(S, 1, x) .* eval_PBSpline.(S, i, x)),
                k * S.h, (k + 1) * S.h, q)
        end
    end
    for j in 2:S.nₕ
        for i in 1:S.nₕ
            M[j, i] = M[j-1, (i-1+S.nₕ-1)%S.nₕ+1]
        end
    end
    return M
end

function stiffnessmatrix(S::PBSpline{T}, q=S.p) where {T}
    K = zeros(T, S.nₕ, S.nₕ)
    for i in 1:S.nₕ
        for k in 0:S.p
            K[1, i] += integrate_gausslegendre(x -> (eval_deriv_PBSpline.(S, 1, x) .* eval_deriv_PBSpline.(S, i, x)),
                k * S.h, (k + 1) * S.h, q)
        end
    end
    for j in 2:S.nₕ
        for i in 1:S.nₕ
            K[j, i] = K[j-1, (i-1+S.nₕ-1)%S.nₕ+1]
        end
    end
    return K
end


struct PoissonSolverPBSplines{DT<:Real} <: PoissonSolver{DT}
    p::Int
    domain::Tuple{DT,DT}

    knots::Vector{DT}
    basis::PBSpline{DT}

    M::Matrix{DT}
    S::Matrix{DT}
    P::Matrix{DT}
    R::Matrix{DT}

    Mfac::LU{DT,Matrix{DT}}
    Sfac::LU{DT,Matrix{DT}}

    ρ::Vector{DT}
    ϕ::Vector{DT}
    rhs::Vector{DT}

    function PoissonSolverPBSplines(domain::Tuple{DT,DT}, order::Int, n::Int) where {DT}
        knots = range(domain[begin], domain[end], n)
        basis = PBSpline(order, n, knots[end] - knots[begin])

        M = massmatrix(basis)
        S = stiffnessmatrix(basis)

        A = ones(n)
        R = A * A' / (A' * A)
        P = Matrix(I, n, n) .- R

        new{DT}(order, domain, knots, basis, M, S, P, R, lu(M), lu(S .+ R),
                zeros(DT, n), zeros(DT, n), zeros(DT, n))
    end
end

Base.length(p::PoissonSolverPBSplines) = length(p.knots)

function solve!(p::PoissonSolverPBSplines, rhs::AbstractVector)
    ldiv!(p.ϕ, p.Sfac, rhs)
    return p
end

function solve!(p::PoissonSolverPBSplines, x::AbstractVector, w::AbstractVector)
    rhs_particles_PBSBasis(x, w, p.basis, p.rhs)
    solve!(p, - p.P * p.rhs)
    ldiv!(p.ρ, p.Mfac, p.rhs)
    return p
end

function solve!(p::PoissonSolverPBSplines, x::AbstractMatrix, w::AbstractMatrix)
    @assert size(x, 1) == 1
    @assert size(w, 1) == 1

    X = reshape(x, size(x, 2))
    W = reshape(w, size(w, 2))

    solve!(p, X, W)
end

function eval_density(p::PoissonSolverPBSplines, x::Real)
    eval_PBSBasis(p.ρ, p.basis, x)
end

function eval_potential(p::PoissonSolverPBSplines, x::Real)
    eval_PBSBasis(p.ϕ, p.basis, x)
end

function eval_field(p::PoissonSolverPBSplines, x::Real)
    eval_deriv_PBSBasis(p.ϕ, p.basis, x)
end
