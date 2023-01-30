
using FFTW
using StatsBase


fftmod(x, domain) = mod(x - domain[begin], domain[end] - domain[begin]) + domain[begin]

struct FFTWBasis{DT, GT <: AbstractVector{DT}}
    domain::Tuple{DT,DT}
    xgrid::GT
    Δx::DT
end

function FFTWBasis(domain, ngrid)
    xgrid = range(domain[begin], domain[end], ngrid + 1)
    Δx = xgrid[2] - xgrid[1]
    FFTWBasis(domain, xgrid, Δx)
end

PeriodicBasisFFTW(domain, n) = FFTWBasis(domain, n)


function get_indices(b::FFTWBasis, x)
    y = fftmod(x, b.domain) - b.domain[begin]

    i1 = floor(Int, y / b.Δx)
    i2 = mod( ceil(Int, y / b.Δx) - 1, length(b)) + 1

    i1 == 0 && (i1 = length(b))
    i2 == 0 && (i2 = length(b))

    i1 == length(b)+1 && (i1 = 1)
    i2 == length(b)+1 && (i2 = 1)

    return (i1, i2)
end

function get_index(b::FFTWBasis, x)
    i1, i2 = get_indices(b, x)
    return (abs(b.xgrid[i1] - x) ≤ abs(b.xgrid[i2] - x) ? i1 : i2)
end


struct FFTWSolution{CT, BT}
    basis::BT
    coefficients::CT
end

function (s::FFTWSolution)(x::Number)
    s.coefficients[get_index(s.basis, x)]
end

PoissonSolution(basis::FFTWBasis, coeffs::AbstractVector) = FFTWSolution(basis, coeffs)

function evalsolution(basis::FFTWBasis, coeffs::AbstractVector, x::Real)
    ϕ = FFTWSolution(basis, coeffs)
    ϕ(x)
end


struct PoissonSolverFFT{DT, BT <: FFTWBasis{DT}} <: PoissonSolver{DT}
    basis::BT
end

Base.length(p::PoissonSolverFFT) = length(p.basis.xgrid) - 1

function solve!(coeffs::AbstractVector, p::PoissonSolverFFT, rhs::AbstractVector)
    ρ̂ = rfft(rhs)
    k² = [(i - 1)^2 for i in eachindex(ρ̂)]
    ϕ̂ = -ρ̂ ./ k²
    ϕ̂[1] = 0
    ϕ̂ ./= ( 2π / (p.basis.domain[end] - p.basis.domain[begin]) )^2
    coeffs .= irfft(ϕ̂, length(rhs))
    return coeffs
end

