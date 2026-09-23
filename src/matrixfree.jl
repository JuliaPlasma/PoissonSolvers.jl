### 1D Laplace operator with periodic bc
function _apply_Δₓ!(y::AbstractVector, x::AbstractVector, h₁)
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())

    @inbounds for i in 1:nx
        i₋ = mod1(i-1, nx)
        i₊ = mod1(i+1, nx)
        y[i] = 1 / (h₁^2) * (x[i₋] - 2*x[i] + x[i₊])
    end
end

function _apply_Δₓ₄!(y::AbstractVector, x::AbstractVector, h₁) # 1D Laplace operator (4th order) with periodic bc
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())
    @inbounds for i in 1:nx
        i₋ = mod1(i-1, nx)
        i₊ = mod1(i+1, nx)
        i₋₋ = mod1(i-2, nx)
        i₊₊ = mod1(i+2, nx)
        y[i] = - 1 / (12 * h₁^2) * (x[i₊₊] - 16*x[i₊] + 30*x[i] - 16*x[i₋] + x[i₋₋])
    end
end

# L = -Δ  + R    
function _apply_Lₓ₄!(y::AbstractVector, x::AbstractVector, h₁) # 1D Laplace operator (4th order) with periodic bc
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())
    Σx = sum(x) / nx
    @inbounds for i in 1:nx
        i₋ = mod1(i-1, nx)
        i₊ = mod1(i+1, nx)
        i₋₋ = mod1(i-2, nx)
        i₊₊ = mod1(i+2, nx)
        y[i] = 1 / (12 * h₁^2) * (x[i₊₊] - 16*x[i₊] + 30*x[i] - 16*x[i₋] + x[i₋₋])
        y[i] += Σx
    end
end

### Constant Nullspace Projection
# if 1 ∈ Δ, then Δϕ = ρ is not well posed but (Δ + R)ϕ = (1 - R)ρ is. 
function _apply_Rₓ!(y::AbstractVector, x::AbstractVector) # Nullspace projection
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())
    Σx = sum(x) / nx
    @inbounds for i in 1:nx
        y[i] = Σx
    end
end

"""
    FiniteDifferenceBasis(domain, ngrid; order = 2)

A uniform periodic grid of `ngrid` cells on `domain`, on which ``-\\Delta`` is the central
difference stencil of the given `order`, 2 or 4.

The grid is the one [`FFTWBasis`](@ref) lays out: `xgrid` holds `ngrid + 1` points, and the last
is not a degree of freedom. This basis selects [`PoissonSolverMatrixFree`](@ref).
"""
struct FiniteDifferenceBasis{DT, GT <: AbstractVector{DT}}
    domain::Tuple{DT, DT}
    xgrid::GT
    Δx::DT
    order::Int
end

function FiniteDifferenceBasis(domain, ngrid; order = 2)
    order ∈ (2, 4) ||
        throw(ArgumentError("the stencil order must be 2 or 4; got order = $(order)"))
    grid = FFTWBasis(domain, ngrid)
    FiniteDifferenceBasis(grid.domain, grid.xgrid, grid.Δx, order)
end

Base.length(b::FiniteDifferenceBasis) = length(b.xgrid) - 1
Base.eltype(::Type{<:FiniteDifferenceBasis{DT}}) where {DT} = DT
Base.eltype(b::FiniteDifferenceBasis) = eltype(typeof(b))
ndofs(b::FiniteDifferenceBasis) = length(b)

# The solution is a grid function on the same grid as the FFT backend's, so it evaluates the
# same way: `FFTWSolution` reads `domain`, `xgrid` and `Δx`, and finds points through these two.
grid(b::FiniteDifferenceBasis) = FFTWBasis(b.domain, b.xgrid, b.Δx)
nearest_indices(b::FiniteDifferenceBasis, x) = nearest_indices(grid(b), x)
nearest_index(b::FiniteDifferenceBasis, x) = nearest_index(grid(b), x)

PoissonSolution(b::FiniteDifferenceBasis, coeffs::AbstractVector) = FFTWSolution(b, coeffs)

@doc raw"""
    PoissonSolverMatrixFree(basis::FiniteDifferenceBasis; reltol, maxiter)

A matrix-free solver for ``-\phi'' = \rho`` on a periodic uniform grid.

The periodic Laplacian annihilates the constants, so ``-\Delta \phi = \rho`` is not well posed.
With ``R`` the projection onto the constants, ``(-\Delta + R) \phi = (1 - R) \rho`` is: the
operator is symmetric positive definite, and its solution is the mean-free one, which is the
solution `PoissonSolverFFT` returns. [`solve!`](@ref) finds it by conjugate gradients, applying
the stencils `_apply_Δₓ!` or `_apply_Δₓ₄!` and `_apply_Rₓ!` and building no matrix.

The iteration stops once the residual norm falls to `reltol` times that of ``(1 - R) \rho``, and
throws an `ErrorException` if `maxiter` iterations do not get it there.

The iteration vectors are fields of the solver, so [`solve!`](@ref) allocates nothing when it is
given a vector of grid values. Given a function, it allocates the vector of samples it takes
first. As for `PoissonSolverFFT`, a solver is therefore not reentrant: give each task its own.
"""
struct PoissonSolverMatrixFree{DT, BT <: FiniteDifferenceBasis{DT}} <: PoissonSolver{DT}
    basis::BT
    reltol::DT
    maxiter::Int
    r::Vector{DT}
    d::Vector{DT}
    Ld::Vector{DT}
    Rx::Vector{DT}

    function PoissonSolverMatrixFree(b::FiniteDifferenceBasis{DT};
            reltol = 16 * eps(DT), maxiter = 4 * length(b)) where {DT}
        r = Vector{DT}(undef, length(b))
        new{DT, typeof(b)}(b, reltol, maxiter, r, similar(r), similar(r), similar(r))
    end
end

PoissonSolver(b::FiniteDifferenceBasis) = PoissonSolverMatrixFree(b)

basis(p::PoissonSolverMatrixFree) = p.basis
Base.length(p::PoissonSolverMatrixFree) = length(p.basis)

gridvalues(p::PoissonSolverMatrixFree, f) = f.(p.basis.xgrid[1:(end - 1)])

# y = (-Δ + R) x
function _apply_L!(y::AbstractVector, x::AbstractVector, p::PoissonSolverMatrixFree)
    b = p.basis
    b.order == 2 ? _apply_Δₓ!(y, x, b.Δx) : _apply_Δₓ₄!(y, x, b.Δx)
    _apply_Rₓ!(p.Rx, x)
    y .= p.Rx .- y
    return y
end

function solve!(ϕ::AbstractVector, p::PoissonSolverMatrixFree, ρ::AbstractVector)
    length(ϕ) == length(ρ) == length(p) || throw(DimensionMismatch(
        "the solver has $(length(p)) degrees of freedom, but the right-hand side has " *
        "$(length(ρ)) and the result $(length(ϕ))"))
    r, d, Ld = p.r, p.d, p.Ld

    # r = (1 - R) ρ is the residual of the initial guess ϕ = 0.
    _apply_Rₓ!(p.Rx, ρ)
    r .= ρ .- p.Rx
    ϕ .= 0
    d .= r
    rr = rr₀ = dot(r, r)
    tolerance = p.reltol^2 * rr₀

    for _ in 1:(p.maxiter)
        rr ≤ tolerance && return ϕ
        _apply_L!(Ld, d, p)
        α = rr / dot(d, Ld)
        ϕ .+= α .* d
        r .-= α .* Ld
        rr, rr₋ = dot(r, r), rr
        d .= r .+ (rr / rr₋) .* d
    end

    rr ≤ tolerance || error("conjugate gradients did not reach a relative residual of " *
          "$(p.reltol) in $(p.maxiter) iterations; it stopped at $(sqrt(rr / rr₀))")
    return ϕ
end

function solve!(ϕ::AbstractVector, p::PoissonSolverMatrixFree, rhs::Base.Callable)
    solve!(ϕ, p, gridvalues(p, rhs))
end

solve(p::PoissonSolverMatrixFree, rhs::AbstractVector) = solve!(similar(rhs), p, rhs)
solve(p::PoissonSolverMatrixFree, rhs::Base.Callable) = solve(p, gridvalues(p, rhs))
