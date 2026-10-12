# Tight solves of every backend: the normwise backward error of the discrete equation and
# conversion of a function right-hand side to the solver's element type.
#
# The file tests three source files — `matrixfree.jl`, `poisson_fft.jl` and `poisson_spline.jl` —
# so it sits in `test/verification/` rather than beside any one of them.
#
# `A`, `b` and `φ` are evaluated in `widen(T)`: Float64 for Float32, and BigFloat for Float64. In
# `T` the residual would carry a stencil's own rounding, and the number bounded here is the
# solver's.

using LinearAlgebra
using PoissonSolvers
using PoissonSolvers: _apply_Δₓ!, _apply_Δₓ₄!, discretise, loadvector
using SimpleSplines
using Test

# The grid of the suite's `Float32` accuracy test; every backend here builds at this size.
const NSOLVE = 256

# A density with every Fourier mode, so that no backend passes on a grid that happens to resolve a
# single one.
density(x) = 4π^2 * exp(sin(2π * x)) * (sin(2π * x) - cos(2π * x)^2)

# The right-hand side of the function `f` in `T`, the caller's own conversion of each sample:
# the grid samples, or the Galerkin load vector. It does not call `discretise`, so it is a
# baseline independent of the hook under test.
function reference(p, ::Type{T}, f) where {T}
    g = x -> T(f(x))
    b = basis(p)
    return b isa AbstractBSplineBasis ? loadvector(p, g) : g.(b.xgrid[1:(end - 1)])
end

# A result of `n` entries in `T`, as a vector or as a contiguous view at an offset. It holds 7,
# away from any solution, so a `solve!` that returns it untouched does not reach the reference.
function container(::Type{T}, n, shape) where {T}
    shape === :vector ? fill(T(7), n) : view(fill(T(7), n + 2), 2:(n + 1))
end

# Eleven backend configurations: FFT, both matrix-free stencil orders, and periodic and
# Dirichlet splines at orders 2 to 5.
function backends(::Type{T}) where {T}
    domain = T.((0.0, 1.0))
    return [("FFT", PoissonSolverFFT(FFTWBasis(domain, NSOLVE)));
            [("matrix-free order $k",
                 PoissonSolverMatrixFree(FiniteDifferenceBasis(domain, NSOLVE; order = k)))
             for k in (2, 4)];
            [("periodic spline order $k",
                 PoissonSolverSpline(PeriodicBasisSpline(domain, k, NSOLVE))) for k in 2:5];
            [("Dirichlet spline order $k",
                 PoissonSolverSpline(DirichletBasisSpline(domain, k, NSOLVE))) for k in 2:5]]
end

# The circulant whose first column is `c`.
function circulant(c)
    n = length(c)
    return [c[mod(i - j, n) + 1] for i in 1:n, j in 1:n]
end

# The first column of the FFT backend's circulant. Its symbol is `s_m = (2π min(m, n − m)/L)²` for
# `m ≥ 1`, and `s_0 = 1`, the k = 0 mode carrying no information. The symbol is real and even,
# so the column is its inverse cosine transform, written out because no FFT of this tree carries
# `widen(T)`.
function fft_column(b, ::Type{W}) where {W}
    n = length(b)
    L = W(b.domain[end]) - W(b.domain[begin])
    s = [m == 0 ? one(W) : (2 * W(π) * min(m, n - m) / L)^2 for m in 0:(n - 1)]
    return [sum(s[m + 1] * cospi(W(2 * m * j) / n) for m in 0:(n - 1)) / n
            for j in 0:(n - 1)]
end

# `−Δ_h + R`, the operator the matrix-free solver inverts: the first column from the package's own
# stencil, circulated, so that `opnorm` sees the matrix the residual does. The periodic stencil is
# shift-invariant, so one column is the whole matrix.
function difference_operator(p::PoissonSolverMatrixFree, ::Type{W}) where {W}
    b = basis(p)
    n = length(b)
    h = W(b.Δx)
    x = zeros(W, n)
    y = similar(x)
    x[1] = one(W)
    b.order == 2 ? _apply_Δₓ!(y, x, h) : _apply_Δₓ₄!(y, x, h)
    return circulant(sum(x) / n .- y)
end

# The discrete operator each backend inverts, dense and in `W`.
function operator(p, ::Type{W}) where {W}
    b = basis(p)
    b isa FFTWBasis && return circulant(fft_column(b, W))
    b isa FiniteDifferenceBasis && return difference_operator(p, W)
    S = stiffness_matrix(p.quadrature)
    if b isa PeriodicBSplineBasis
        # `SimpleSplines.CirculantMass` inverts the circulant of the first stiffness column,
        # not the assembled matrix, which can deviate from circulant through rounding.
        # Adding the mean projection makes the constant mode nonsingular.
        return circulant(W.(Matrix(S)[:, 1])) .+ one(W) / length(p)
    end
    return W.(Matrix(S))
end

# The kernel of each periodic operator holds the constants, so its right-hand side is mean-free.
periodic_rhs(p) = !(basis(p) isa RecombinedBSplineBasis)

# The normwise backward error of Rigal and Gaches: `(A + ΔA)φ = b + Δb` with `‖ΔA‖∞ ≤ η‖A‖∞` and
# `‖Δb‖∞ ≤ η‖b‖∞` (Higham, *Accuracy and Stability of Numerical Algorithms*, Thm 7.1).
function backward_error(A, b, φ)
    norm(b - A * φ, Inf) / (opnorm(A, Inf) * norm(φ, Inf) + norm(b, Inf))
end

@testset "the normwise backward error of the discrete equation" begin
    # 16 eps(T) is the matrix-free backend's default `reltol`, and it is written as the literal so
    # that a looser default fails here.
    for T in (Float32, Float64)
        W = widen(T)
        for (name, p) in backends(T)
            A = operator(p, W)
            periodic = periodic_rhs(p)
            n = length(p)
            @testset "$name, $T" begin
                for scale in (1.0, 3e-8, 3e8), shape in (:vector, :view)
                    # The scaled right-hand side reaches the solver the same way through a
                    # contiguous view at an offset; the backward error is scale-invariant.
                    ρ = reference(p, T, x -> scale * density(x))
                    rhs = shape === :vector ? ρ : copyto!(container(T, n, shape), ρ)
                    φ = solve!(similar(ρ), p, rhs)
                    b = W.(ρ)
                    periodic && (b .-= sum(b) / n)
                    @test backward_error(A, b, W.(φ)) ≤ 16 * eps(T)
                end
            end
        end
    end
end

# A right-hand side whose values are never the solver's element type: `Float64` values at a
# `Float32` point, through the `Float64` product `2π * x`, and `Float32` values at a `Float64` one.
other(x::Float32) = exp(sin(2π * x))
other(x::Float64) = Float32(exp(sin(2π * x)))

# The least of five warmed readings of the bytes `g()` allocates.
function allocation(g::G) where {G}
    g()
    bytes = typemax(Int)
    for _ in 1:5
        bytes = min(bytes, @allocated g())
    end
    return bytes
end

# The slack of an allocation comparison. One `@allocated` reading that reaches `exp` or `inv`
# varies by 16–64 B between calls on Windows; a copy of the right-hand side, which the comparison
# has to catch, is at least `NSOLVE * sizeof(Float32)`, 1 KiB.
const SLACK = 128

# Check the element type, sample conversion, preallocated result and allocations for one backend.
# `T` is a type parameter, so `@inferred` sees a fixed type.
function conversion_probe(p, f, ::Type{T}) where {T}
    n = length(p)
    ρ = reference(p, T, f)
    φ = @inferred solve(p, f)
    @test eltype(φ) === T
    @test φ == solve(p, ρ)

    # Reducing `f` costs no more than the caller's own conversion, whatever `f` returns.
    @test allocation(() -> discretise(p, f)) ≤ allocation(() -> reference(p, T, f)) + SLACK

    for shape in (:vector, :view)
        out = container(T, n, shape)
        @test solve!(out, p, f) === out

        # Against a result of the same layout: a view and a vector can reach different BLAS
        # kernels, which agree only to rounding.
        @test out == solve!(container(T, n, shape), p, ρ)

        # The whole call, not only the reduction above: a function path that copies the
        # right-hand side on its way into `out` costs one buffer more than the caller's own
        # conversion does.
        @test allocation(() -> solve!(out, p, f)) ≤
              allocation(() -> solve!(out, p, reference(p, T, f))) + SLACK
    end
    return nothing
end

@testset "the element type of a function right-hand side" begin
    # Both directions, so that a solver converting to a fixed `Float32` fails the `Float64` half.
    for T in (Float32, Float64), (name, p) in backends(T)

        @testset "$name, $T" begin
            conversion_probe(p, other, T)
        end
    end
end
