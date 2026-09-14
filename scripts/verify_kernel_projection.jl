# The null space of the periodic stiffness matrix, and what it costs to remove it.
#
# Run in a fresh process:
#
#     julia --startup-file=no --project=. scripts/verify_kernel_projection.jl
#
# Constants lie in the kernel of the periodic stiffness matrix, so `-φ'' = ρ` fixes φ only up to
# one. Two cures give the same answer. The rank-one shift by the mean projector, `S + 𝟙𝟙ᵀ/N`,
# is invertible and is what a textbook reaches for; the deflation this package uses leaves the
# matrix alone and solves in the complement of the kernel instead.
#
# The claims established here:
#
#   1. The kernel is the constants, on every mesh, and the transform of a uniform assembly puts
#      that kernel in one Fourier coefficient.
#   2. The shift changes that coefficient from zero to one and changes nothing else. The O(N²)
#      fill-in it costs therefore buys exactly one number.
#   3. The deflated solve reproduces the shifted one to round-off, on both representations, and
#      returns the mean-free solution exactly rather than to within the shift's own rounding.
#   4. Storage and construction stay O(N) rather than O(N²).
#   5. The shift is an absolute constant while the spectrum is not, so on a domain that is not
#      of order one it also costs two orders of condition number. The deflation has no scale to
#      get wrong.

using FFTW
using LinearAlgebra
using PoissonSolvers
using Printf
using SimpleSplines
using SparseArrays
using Test

const DOMAIN = (0.0, 2π)
const ORDER = 5                        # k = p + 1, so the degree is 4

basis(meshtype, n) = PeriodicBSplineBasis(meshtype(n, DOMAIN[2] - DOMAIN[1]), ORDER - 1)

"""The stiffness matrix of a periodic basis of `n` cells on `meshtype`."""
stiffness(meshtype, n) = stiffness_matrix(SplineQuadrature(basis(meshtype, n)))

header(s) = (println(); println(s); println("-"^length(s)))

# ---------------------------------------------------------------------------------------------
function kernel_is_the_constants()
    header("1. the kernel is the constants, and a uniform assembly puts it in one coefficient")

    for meshtype in (UniformMesh, GradedMesh, RandomMesh), n in (32, 128, 512)

        S = stiffness(meshtype, n)
        residual = norm(S * ones(n), Inf) / norm(S, Inf)
        @test residual < 1e-12

        if meshtype === UniformMesh
            ĉ = rfft(Vector(S[:, 1]))
            @printf("%-11s n=%3d  ‖S𝟙‖/‖S‖=%.2e  |ĉ₀|/max|ĉ|=%.2e  next smallest=%.2e\n",
                nameof(meshtype), n, residual,
                abs(ĉ[1]) / maximum(abs, ĉ),
                minimum(abs, ĉ[2:end]) / maximum(abs, ĉ))
            @test abs(ĉ[1]) / maximum(abs, ĉ) < 1e-13
            @test minimum(abs, ĉ[2:end]) / maximum(abs, ĉ) > 1e-6
        else
            @printf("%-11s n=%3d  ‖S𝟙‖/‖S‖=%.2e  (not circulant, so no transform)\n",
                nameof(meshtype), n, residual)
        end
    end
end

# ---------------------------------------------------------------------------------------------
function the_shift_changes_one_coefficient()
    header("2. the rank-one shift changes that one coefficient and nothing else")

    for n in (32, 128, 512)
        S = stiffness(UniformMesh, n)
        c = Vector(S[:, 1])
        e₀ = [k == 1 ? one(ComplexF64) : zero(ComplexF64) for k in 1:(n ÷ 2 + 1)]

        # 𝟙𝟙ᵀ/N is circulant too, with first column 𝟙/N, and its transform is (1, 0, 0, …)
        defect = norm(rfft(c .+ inv(n)) .- rfft(c) .- e₀, Inf)
        @printf("n=%3d  ‖ℱ(c + 𝟙/N) - ℱ(c) - e₀‖∞ = %.2e   stored entries %6d → %7d\n",
            n, defect, nnz(S), nnz(S .+ inv(n)))
        @test defect < 1e-14
    end
end

# ---------------------------------------------------------------------------------------------
function the_deflation_reproduces_the_shift()
    header("3. the deflated solve is the shifted solve, on both representations")

    for meshtype in (UniformMesh, GradedMesh), n in (32, 127)

        b = basis(meshtype, n)
        S = stiffness(meshtype, n)
        op = mass_operator(S, b; kernel = :project)
        shifted = mass_operator(Matrix(S) .+ inv(n), b)

        worst, gauge = 0.0, 0.0
        for _ in 1:5
            rhs = randn(n)
            meanfree = rhs .- sum(rhs) / n
            φ = op \ rhs

            worst = max(worst, maximum(abs, φ .- shifted \ meanfree))
            gauge = max(gauge, abs(sum(φ) / n))
            @test S * φ ≈ meanfree atol = 1e-10
            @test op \ meanfree ≈ φ atol = 1e-12     # the mean is dropped either way
        end

        @printf("%-11s n=%3d  %-14s ‖deflated - shifted‖∞=%.2e  |mean|=%.2e\n",
            nameof(meshtype), n, nameof(typeof(op)), worst, gauge)
        @test worst < 1e-10
        @test gauge < 1e-14
    end
end

# ---------------------------------------------------------------------------------------------
function the_solver_stays_linear()
    header("4. the solver stays O(N)")

    for n in (256, 1024, 2048)
        b = basis(UniformMesh, n)
        S = stiffness(UniformMesh, n)

        # `S .+ inv(n)` is what the shift was written as, and it stays a SparseMatrixCSC —
        # structurally full, so it carries a row index beside every one of the N² entries.
        deflated = Base.summarysize(PoissonSolverSpline(b))
        shift = Base.summarysize(S .+ inv(n))
        @printf("n=%4d  whole solver %8d B   the shifted matrix alone %9d B  (%.0f×)\n",
            n, deflated, shift, shift / deflated)
        @test deflated < shift
    end
end

# ---------------------------------------------------------------------------------------------
function the_shift_picks_a_scale()
    header("5. the shift is an absolute constant, the spectrum is not")

    n = 64
    for length_ in (2π * 1e-3, 2π, 2π * 1e3)
        b = PeriodicBSplineBasis(UniformMesh(n, length_), ORDER - 1)
        S = Matrix(stiffness_matrix(SplineQuadrature(b)))
        λ = sort(eigvals(Symmetric(S)))
        μ = sort(eigvals(Symmetric(S .+ inv(n))))

        @printf("L=%9.3e  κ on the mean-free subspace=%.3e  κ(S + 𝟙𝟙ᵀ/N)=%.3e  (%.0f×)\n",
            length_, λ[end] / λ[2], μ[end] / μ[1], (μ[end] / μ[1]) / (λ[end] / λ[2]))
    end
    println("\nThe deflation works on the mean-free subspace itself, so it has no such scale.")
end

function report()
    kernel_is_the_constants()
    the_shift_changes_one_coefficient()
    the_deflation_reproduces_the_shift()
    the_solver_stays_linear()
    the_shift_picks_a_scale()
    println()
end

report()
