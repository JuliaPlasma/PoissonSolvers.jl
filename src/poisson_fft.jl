using FFTW
using LinearAlgebra

fftmod(x, domain) = mod(x - domain[begin], domain[end] - domain[begin]) + domain[begin]

"""
    FFTWBasis(domain, ngrid)

A uniform periodic grid of `ngrid` cells on `domain`.

`xgrid` holds `ngrid + 1` points, both endpoints included; the last coincides with the first
under periodicity and is not a degree of freedom, so `length(b) == ngrid`.
"""
struct FFTWBasis{DT, GT <: AbstractVector{DT}}
    domain::Tuple{DT, DT}
    xgrid::GT
    Δx::DT
end

function FFTWBasis(domain, ngrid)
    xgrid = range(domain[begin], domain[end], ngrid + 1)
    FFTWBasis((domain[begin], domain[end]), xgrid, step(xgrid))
end

Base.length(b::FFTWBasis) = length(b.xgrid) - 1
ndofs(b::FFTWBasis) = length(b)

function nearest_indices(b::FFTWBasis, x)
    i1 = floor(Int, (x - b.domain[begin]) / b.Δx) + 1
    i2 = i1 + 1
    return (mod(i1 - 1, length(b)) + 1, mod(i2 - 1, length(b)) + 1)
end

function nearest_index(b::FFTWBasis, x)
    i1, i2 = nearest_indices(b, x)
    i = (abs(b.xgrid[i1] - fftmod(x, b.domain)) ≤ abs(b.xgrid[i2] - fftmod(x, b.domain)) ?
         i1 : i2)
    return i
end

"""
    FFTWSolution(basis, coefficients)

The solution on an [`FFTWBasis`](@ref), evaluated at the nearest grid point.
"""
struct FFTWSolution{CT, BT}
    basis::BT
    coefficients::CT
end

basis(s::FFTWSolution) = s.basis
coefficients(s::FFTWSolution) = s.coefficients

(s::FFTWSolution)(x::Number) = s.coefficients[nearest_index(s.basis, x)]

function (s::FFTWSolution)(x::Number, d::Integer)
    d == 0 && return s(x)
    d == 1 || throw(ArgumentError(
        "a grid solution carries only a first derivative, a one-sided difference; got d = $(d)"))
    i1, i2 = nearest_indices(s.basis, x)
    return (s.coefficients[i2] - s.coefficients[i1]) / s.basis.Δx
end

"""
    FFTWDerivative

A derivative of an [`FFTWSolution`](@ref), as returned by `derivative`. Callable, so that
`derivative(s).(v)` broadcasts the way `s.(v)` does, matching the spline side.
"""
struct FFTWDerivative{ST <: FFTWSolution, OT}
    solution::ST
    d::OT
end

(ds::FFTWDerivative)(x::Number) = ds.solution(x, ds.d)

derivative(s::FFTWSolution, d = 1) = FFTWDerivative(s, d)

PoissonSolution(b::FFTWBasis, coeffs::AbstractVector) = FFTWSolution(b, coeffs)

"""
    PoissonSolverFFT(basis)

A spectral solver for ``-\\phi'' = \\rho`` on a periodic uniform grid.

The transforms and the inverse Laplacian symbol are built once, so [`solve!`](@ref) allocates
nothing when it is given a vector of grid values. Given a function, it allocates the vector of
samples it takes first.
"""
struct PoissonSolverFFT{DT, BT <: FFTWBasis{DT}, PT, IT} <: PoissonSolver{DT}
    basis::BT
    plan::PT
    iplan::IT
    k⁻²::Vector{DT}
    ρ̂::Vector{Complex{DT}}

    function PoissonSolverFFT(b::FFTWBasis{DT}) where {DT}
        n = length(b)
        ρ̂ = Vector{Complex{DT}}(undef, n ÷ 2 + 1)

        # UNALIGNED, so that the plans accept any strided argument — a view into a larger array
        # in particular, whose alignment an aligned plan rejects at run time. ESTIMATE has to be
        # given alongside it: UNALIGNED alone replaces the flags rather than adding to them, and
        # FFTW's default rigor then measures, which overwrites the array being planned for.
        flags = FFTW.ESTIMATE | FFTW.UNALIGNED
        plan = plan_rfft(Vector{DT}(undef, n); flags)
        iplan = plan_irfft(similar(ρ̂), n; flags)

        # -φ'' = ρ is k² φ̂ = ρ̂ with k = 2πm/L. The m = 0 mode is the constant, which the
        # equation does not determine at all; a zero factor there picks the mean-free solution
        # rather than dividing by zero and overwriting the result afterwards.
        L = b.domain[end] - b.domain[begin]
        k⁻² = DT[m == 0 ? zero(DT) : inv((2π * m / L)^2) for m in 0:(length(ρ̂) - 1)]

        new{DT, typeof(b), typeof(plan), typeof(iplan)}(b, plan, iplan, k⁻², ρ̂)
    end
end

PoissonSolver(b::FFTWBasis) = PoissonSolverFFT(b)

basis(p::PoissonSolverFFT) = p.basis
Base.length(p::PoissonSolverFFT) = length(p.basis)

gridvalues(p::PoissonSolverFFT, f) = f.(p.basis.xgrid[1:(end - 1)])

function solve!(coeffs::AbstractVector, p::PoissonSolverFFT, rhs::AbstractVector)
    mul!(p.ρ̂, p.plan, rhs)
    p.ρ̂ .*= p.k⁻²
    mul!(coeffs, p.iplan, p.ρ̂)
    return coeffs
end

function solve!(coeffs::AbstractVector, p::PoissonSolverFFT, rhs::Base.Callable)
    solve!(coeffs, p, gridvalues(p, rhs))
end

solve(p::PoissonSolverFFT, rhs::AbstractVector) = solve!(similar(rhs), p, rhs)
solve(p::PoissonSolverFFT, rhs::Base.Callable) = solve(p, gridvalues(p, rhs))
