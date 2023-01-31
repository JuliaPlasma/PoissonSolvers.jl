
using BSplineKit.Splines: PeriodicVector
using FFTW


fftmod(x, domain) = mod(x - domain[begin], domain[end] - domain[begin]) + domain[begin]


struct FFTWBasis{DT, GT <: AbstractVector{DT}}
    domain::Tuple{DT,DT}
    xgrid::PeriodicVector{DT,GT}
    Δx::DT
end

function FFTWBasis(domain, ngrid)
    xgrid = PeriodicVector(range(domain[begin], domain[end], ngrid + 1))
    Δx = xgrid[2] - xgrid[1]
    FFTWBasis(domain, xgrid, Δx)
end

PeriodicBasisFFTW(domain, n) = FFTWBasis(domain, n)

Base.length(basis::FFTWBasis) = length(basis.xgrid) - 1


function nearest_indices(b::FFTWBasis, x)
    i1 = floor(Int, (x - b.domain[begin]) / b.Δx) + 1
    i2 = i1 + 1
    return (mod(i1-1, length(b)) + 1, mod(i2-1, length(b)) + 1)
end

function nearest_index(b::FFTWBasis, x)
    i1, i2 = nearest_indices(b, x)
    i = (abs(b.xgrid[i1] - fftmod(x, b.domain)) ≤ abs(b.xgrid[i2] - fftmod(x, b.domain)) ? i1 : i2)
    return i
end


struct FFTWSolution{CT, BT}
    basis::BT
    coefficients::CT
end

function (s::FFTWSolution)(x::Number)
    s.coefficients[nearest_index(s.basis, x)]
end

PoissonSolution(basis::FFTWBasis, coeffs::AbstractVector) = FFTWSolution(basis, coeffs)

function evalsolution(basis::FFTWBasis, coeffs::AbstractVector, x::Real)
    ϕ = FFTWSolution(basis, coeffs)
    ϕ(x)
end


struct FFTWDerivative{N, ST <: FFTWSolution}
    solution::ST
    function FFTWDerivative{N}(solution::ST) where {N,ST}
        new{N,ST}(solution)
    end
end

Base.:*(::Derivative{0}, s::FFTWSolution) = s
Base.:*(::Derivative{N}, s::FFTWSolution) where {N} = FFTWDerivative{N}(s)

Base.diff(s::FFTWSolution, op = Derivative(1)) = op * s

function (d::FFTWDerivative{1})(x::Number)
    i1, i2 = nearest_indices(d.solution.basis, x)
    return (d.solution.coefficients[i2] - d.solution.coefficients[i1]) / d.solution.basis.Δx
end


struct PoissonSolverFFT{DT, BT <: FFTWBasis{DT}} <: PoissonSolver{DT}
    basis::BT
end

Base.length(p::PoissonSolverFFT) = length(p.basis)

PoissonSolver(basis::FFTWBasis) = PoissonSolverFFT(basis)

function solve!(coeffs::AbstractVector, p::PoissonSolverFFT, rhs::AbstractVector)
    ρ̂ = rfft(rhs)
    k² = [(i - 1)^2 for i in eachindex(ρ̂)]
    ϕ̂ = ρ̂ ./ k²
    ϕ̂[1] = 0
    ϕ̂ ./= ( 2π / (p.basis.domain[end] - p.basis.domain[begin]) )^2
    coeffs .= irfft(ϕ̂, length(rhs))
    return coeffs
end

function solve!(coeffs::AbstractVector, p::PoissonSolverFFT, rhs::Base.Callable)
    solve!(coeffs, p, rhs.(p.basis.xgrid[1:end-1]))
end

function solve(p::PoissonSolverFFT, rhs::AbstractVector)
    solve!(zero(rhs), p, rhs)
end

function solve(p::PoissonSolverFFT, rhs::Base.Callable)
    solve(p, rhs.(p.basis.xgrid[1:end-1]))
end
