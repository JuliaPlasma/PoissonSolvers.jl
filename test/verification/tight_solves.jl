# The `Float32` solves of every backend on one grid: the normwise backward error of the discrete
# equation (T2.4's *Done when*, §5V.P41 check 1), and the element type a function right-hand side
# is reduced in (check 3).
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

# The grid every backend is checked on. It is the size the suite's `Float32` accuracy is measured
# at, and the largest at which a `Float32` periodic spline basis still builds.
const NSOLVE = 256

# A density with every Fourier mode, so that no backend passes on a grid that happens to resolve a
# single one.
density(x) = 4π^2 * exp(sin(2π * x)) * (sin(2π * x) - cos(2π * x)^2)

# `f`'s samples in `T`, the caller's own conversion: the baseline the solver's conversion is
# measured against. The closure carries both as type parameters, so it is a function barrier of
# its own.
in_element(f::F, ::Type{T}) where {F, T} = x -> T(f(x))

# The right-hand side of the function `f` in `T`: the grid samples, or the Galerkin load vector.
function reference(p, ::Type{T}, f) where {T}
    basis(p) isa AbstractBSplineBasis ? loadvector(p, in_element(f, T)) :
    T.(f.(basis(p).xgrid[1:(end - 1)]))
end

# The eleven backend configurations check 1 names: the FFT backend, the matrix-free backend at
# both stencil orders, and the periodic and Dirichlet spline backends at orders 2 to 5.
function backends(::Type{T}) where {T}
    domain = T.((0.0, 1.0))
    return (("FFT", PoissonSolverFFT(FFTWBasis(domain, NSOLVE))),
        ("matrix-free order 2",
            PoissonSolverMatrixFree(FiniteDifferenceBasis(domain, NSOLVE; order = 2))),
        ("matrix-free order 4",
            PoissonSolverMatrixFree(FiniteDifferenceBasis(domain, NSOLVE; order = 4))),
        ("periodic spline order 2",
            PoissonSolverSpline(PeriodicBasisSpline(domain, 2, NSOLVE))),
        ("periodic spline order 3",
            PoissonSolverSpline(PeriodicBasisSpline(domain, 3, NSOLVE))),
        ("periodic spline order 4",
            PoissonSolverSpline(PeriodicBasisSpline(domain, 4, NSOLVE))),
        ("periodic spline order 5",
            PoissonSolverSpline(PeriodicBasisSpline(domain, 5, NSOLVE))),
        ("Dirichlet spline order 2",
            PoissonSolverSpline(DirichletBasisSpline(domain, 2, NSOLVE))),
        ("Dirichlet spline order 3",
            PoissonSolverSpline(DirichletBasisSpline(domain, 3, NSOLVE))),
        ("Dirichlet spline order 4",
            PoissonSolverSpline(DirichletBasisSpline(domain, 4, NSOLVE))),
        ("Dirichlet spline order 5",
            PoissonSolverSpline(DirichletBasisSpline(domain, 5, NSOLVE))))
end

# The circulant whose first column is `c`.
function circulant(c)
    n = length(c)
    return [c[mod(i - j, n) + 1] for i in 1:n, j in 1:n]
end

# The first column of the FFT backend's circulant. Its symbol is `s_m = (2π min(m, n − m)/L)²` for
# `m ≥ 1`, and `s_0 = 1`, the k = 0 mode carrying no information. The half-spectrum is real and
# even, so the column is the even inverse transform written out: no FFT of this tree carries
# `widen(T)`.
function fft_column(b, ::Type{W}) where {W}
    n = length(b)
    L = W(b.domain[end]) - W(b.domain[begin])
    s = W[(m == 0 ? one(W) : (2 * W(π) * min(m, n - m) / L)^2) for m in 0:(n - 1)]
    c = zeros(W, n)
    for j in 0:(n - 1)
        acc = s[1] + (isodd(j) ? -s[n ÷ 2 + 1] : s[n ÷ 2 + 1])
        for m in 1:(n ÷ 2 - 1)
            acc += 2 * s[m + 1] * cos(2 * W(π) * m * j / n)
        end
        c[j + 1] = acc / n
    end
    return c
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
        # What `SimpleSplines.CirculantMass` inverts is the circulant of `M[:, 1]`, not the
        # assembled stiffness matrix, which deviates from circulant by up to 457 eps; the mean
        # projection is what makes that circulant nonsingular.
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
            @testset "$name" begin
                for scale in (1.0, 3e-8, 3e8), shape in (:vector, :view)

                    ρ = reference(p, T, x -> scale * density(x))
                    # The scaled right-hand side reaches the solver the same way through a
                    # contiguous view at an offset; the backward error is scale-invariant.
                    rhs = if shape === :view
                        buffer = zeros(T, n + 2)
                        buffer[2:(n + 1)] .= ρ
                        view(buffer, 2:(n + 1))
                    else
                        ρ
                    end
                    φ = solve!(similar(ρ), p, rhs)
                    b = W.(ρ)
                    periodic && (b .-= sum(b) / n)
                    @test backward_error(A, b, W.(φ)) ≤ 16 * eps(T)
                end
            end
        end
    end
end

# `f` returns `Float64` for a `Float32` argument, through the `Float64` product `2π * x`.
float64_samples(x) = exp(sin(2π * x))

# `f` returns `Float32` for a `Float64` argument.
float32_samples(x) = Float32(exp(sin(2π * x)))

# The function whose return type is never the solver's element type, for each of the two.
other_element(::Type{Float32}) = float64_samples
other_element(::Type{Float64}) = float32_samples

# The bytes the samples of a function right-hand side cost, through a function barrier. The result
# is preallocated, so only the samples are at stake — and both sides are then one call that returns
# a vector, so that no frame of the `solve!` path is counted on one side only.
function samples_allocation(p, f)
    discretise(p, f)
    return @allocated discretise(p, f)
end

# The unfixed tree's own conversion, which is the baseline the clause names: the caller converts
# the samples in `T` themselves. It must not go through the tree's hook — under a conversion
# written as `T.(f.(x))` the hook moves with the mutant, and the comparison then measures the
# mutant against itself.
function typed_samples(p, f, ::Type{T}) where {T}
    b = basis(p)
    return b isa AbstractBSplineBasis ? loadvector(p, in_element(f, T)) :
           in_element(f, T).(b.xgrid[1:(end - 1)])
end

function baseline_allocation(p, f, ::Type{T}) where {T}
    typed_samples(p, f, T)
    return @allocated typed_samples(p, f, T)
end

# The bytes the whole preallocated function path costs: `solve!` reduces `f` and moves the result
# into the caller's own `out`. `out` is one of the caller's buffers on both sides of the
# comparison below, so the movement cancels and what is left is the right-hand side. The least of
# five warmed readings, because an `@allocated` that reaches `exp` or `inv` differs by 16–48 B
# between calls on Windows, and CI runs `windows-latest`.
function complete_allocation(out::O, p::P, f::F) where {O, P, F}
    solve!(out, p, f)
    bytes = typemax(Int)
    for _ in 1:5
        bytes = min(bytes, @allocated solve!(out, p, f))
    end
    return bytes
end

# The same path on the unfixed tree — the caller makes the samples in `T` themselves and `solve!`
# moves them into `out` — so that the two readings differ by the reduction alone. It calls neither
# `discretise` nor `reference`, so a hook written as `T.(f.(x))` cannot move with the expression it
# is measured against.
function typed_allocation(out::O, p::P, f::F, ::Type{T}) where {O, P, F, T}
    solve!(out, p, typed_samples(p, f, T))
    bytes = typemax(Int)
    for _ in 1:5
        bytes = min(bytes, @allocated solve!(out, p, typed_samples(p, f, T)))
    end
    return bytes
end

# Check 3 for one backend and one `f`: the element type, the reduction to it, the movement of a
# preallocated result, and the allocation of the samples. `T` is a type parameter, so `@inferred`
# sees a fixed one. Each run covers the `Vector` the solver returns and the unit-stride view of
# check 1, whose result the function path must reach too.
function conversion_probe(p, f, ::Type{T}) where {T}
    n = length(p)
    φ = @inferred solve(p, f)
    @test eltype(φ) === T
    @test φ == solve(p, reference(p, T, f))

    # Reducing `f` costs no more than the caller's own conversion, whatever `f` returns.
    @test samples_allocation(p, f) ≤ baseline_allocation(p, f, T)

    for shape in (:vector, :view)
        # The buffer starts away from the solution values, so a `solve!` that returns it untouched
        # does not reach the reference.
        buffer = fill(T(7), shape === :vector ? n : n + 2)
        out = shape === :vector ? buffer : view(buffer, 2:(n + 1))

        @test solve!(out, p, f) === out
        @test out == solve(p, reference(p, T, f))

        # The whole call, not only the reduction above: a function path that copies the
        # right-hand side on its way into `out` costs one buffer more than the caller's own
        # conversion does.
        @test complete_allocation(out, p, f) ≤ typed_allocation(out, p, f, T)
    end
    return nothing
end

@testset "the element type of a function right-hand side" begin
    # Both directions, so that a solver converting to a fixed `Float32` fails the `Float64` half:
    # the right-hand side of a `Float32` solver returns `Float64`, and that of a `Float64` solver
    # returns `Float32`.
    for T in (Float32, Float64)
        f = other_element(T)
        for (name, p) in backends(T)
            @testset "$name" begin
                conversion_probe(p, f, T)
            end
        end
    end
end
