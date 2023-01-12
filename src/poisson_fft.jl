
using FFTW
using StatsBase

fftmod(x, domain) = mod(x - domain[begin], domain[end] - domain[begin]) + domain[begin]

struct PoissonSolverFFT{DT <: Real} <: PoissonSolver{DT}
    domain::Tuple{DT,DT}
    xgrid::Vector{DT}
    Δx::DT

    ρ::Vector{DT}
    ϕ::Vector{DT}

    function PoissonSolverFFT(domain::Tuple{DT,DT}, n::Int) where {DT}
        xgrid = range(domain[begin], domain[end], n + 1)
        Δx = xgrid[2] - xgrid[1]
        new{DT}(domain, xgrid, Δx, zeros(DT, n), zeros(DT, n))
    end
end

Base.length(p::PoissonSolverFFT) = length(p.xgrid) - 1

function solve!(p::PoissonSolverFFT, rhs::AbstractVector)
    p.ρ .= rhs
    ρ̂ = rfft(p.ρ)
    k² = [(i - 1)^2 for i in eachindex(ρ̂)]
    ϕ̂ = -ρ̂ ./ k²
    ϕ̂[1] = 0
    ϕ̂ ./= ( 2π / (p.domain[end] - p.domain[begin]) )^2
    p.ϕ .= irfft(ϕ̂, length(p.ρ))
    return p
end

function solve!(p::PoissonSolverFFT, points::AbstractVector, weights::AbstractVector{DT}) where {DT}
    h = fit(Histogram, map(x -> fftmod(x, p.domain), points), Weights(weights), p.xgrid)
    d = normalize(h; mode = :density)
    solve!(p, DT.(d.weights))
end


function get_indices(p::PoissonSolverFFT, x)
    y = fftmod(x, p.domain) - p.domain[begin]

    i1 = floor(Int, y / p.Δx)
    i2 = mod( ceil(Int, y / p.Δx) - 1, length(p)) + 1

    i1 == 0 && (i1 = length(p))
    i2 == 0 && (i2 = length(p))

    i1 == length(p)+1 && (i1 = 1)
    i2 == length(p)+1 && (i2 = 1)

    return (i1, i2)
end

function get_index(p::PoissonSolverFFT, x)
    i1, i2 = get_indices(p, x)
    return (abs(p.xgrid[i1] - x) ≤ abs(p.xgrid[i2] - x) ? i1 : i2)
end


function eval_density(p::PoissonSolverFFT{DT}, x::DT) where {DT}
    p.ρ[get_index(p, x)]
end

function eval_potential(p::PoissonSolverFFT{DT}, x::DT) where {DT}
    p.ϕ[get_index(p, x)]
end

function eval_field(p::PoissonSolverFFT{DT}, x::DT) where {DT}
    i1, i2 = get_indices(p, x)
    return - (p.ϕ[i2] - p.ϕ[i1]) / p.Δx
end
