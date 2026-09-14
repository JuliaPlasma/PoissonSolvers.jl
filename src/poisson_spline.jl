using LinearAlgebra
using SimpleSplines

@doc raw"""
    PeriodicBasisSpline(domain, order, ncells)

A periodic B-spline basis of the given `order` on `ncells` uniform cells of `domain`.

`order` is ``k = p + 1``, so `order = 2` is the piecewise-linear basis. Note that
`SimpleSplines` constructors take the degree ``p`` instead.
"""
function PeriodicBasisSpline(domain, order, ncells)
    PeriodicBSplineBasis(UniformMesh(ncells, domain), order - 1)
end

@doc raw"""
    DirichletBasisSpline(domain, order, ncells)

A B-spline basis of the given `order` on `ncells` uniform cells of `domain`, recombined so that
every basis function vanishes at both ends. It spans the solutions of ``-\phi'' = \rho`` with
homogeneous Dirichlet boundaries, and has two degrees of freedom fewer than the clamped basis.

`order` is ``k = p + 1``, as for [`PeriodicBasisSpline`](@ref).
"""
function DirichletBasisSpline(domain, order, ncells)
    BSplineBasis(UniformMesh(ncells, domain), order - 1, Dirichlet())
end

# Constants lie in the kernel of the periodic stiffness matrix, so it is singular. Shifting it by
# the rank-one mean projector ``R = 𝟙𝟙ᵀ/n`` makes it invertible, and taking the mean out of the
# right-hand side makes the solution of the shifted system the one we want: ``S𝟙 = 0`` is what
# keeps the constant mode and the rest from mixing. A Dirichlet basis has no constant mode, so
# both operations are the identity there.
regularise(S, ::PeriodicBSplineBasis) = S .+ inv(size(S, 1))
regularise(S, ::AbstractBSplineBasis) = S

meanfree!(y, x, ::PeriodicBSplineBasis) = y .= x .- sum(x) / length(x)
meanfree!(y, x, ::AbstractBSplineBasis) = y .= x

# The kernel of the stiffness matrix holds exactly the constants the basis represents, which is
# what `polynomial_reproduction ≥ 0` reports. The periodic basis is the one such case treated
# here, by the shift above; a clamped basis and a Neumann recombination are singular with no
# treatment. Saying so here is what makes the reason reach the caller: the factorisation further
# down rejects the same matrix, but reports it against the mass matrix and the quadrature order,
# neither of which is what went wrong.
function checkbasis(b::AbstractBSplineBasis)
    polynomial_reproduction(b) ≥ 0 && throw(ArgumentError(
        "the stiffness matrix of this basis is singular: the basis represents the constants, " *
        "which -φ'' = ρ leaves undetermined. Use PeriodicBasisSpline for a periodic problem, " *
        "or DirichletBasisSpline for one with homogeneous Dirichlet boundaries."))
    return nothing
end

checkbasis(::PeriodicBSplineBasis) = nothing

@doc raw"""
    PoissonSolverSpline(basis)

A B-spline Galerkin solver for ``-\phi'' = \rho`` on `basis`.

The stiffness matrix is factorised once, through the representation `SimpleSplines` chooses for
the basis: an FFT for a periodic uniform basis, a banded Cholesky for a Dirichlet one. Both make
[`solve!`](@ref) allocation-free when it is given a coefficient vector. Given a function, it
allocates the load vector it samples first — see [`loadvector`](@ref).

The basis must be periodic or Dirichlet-recombined. Any other basis represents the constants,
which leaves the stiffness matrix singular; the constructor rejects it.
"""
struct PoissonSolverSpline{DT, QT <: SplineQuadrature{DT},
    MT <: MassOperator{DT}} <: PoissonSolver{DT}
    quadrature::QT
    stiffness::MT

    function PoissonSolverSpline(b::AbstractBSplineBasis{DT}) where {DT}
        checkbasis(b)
        q = SplineQuadrature(b)
        S = mass_operator(regularise(stiffness_matrix(q), b), b)
        new{DT, typeof(q), typeof(S)}(q, S)
    end
end

PoissonSolver(b::AbstractBSplineBasis) = PoissonSolverSpline(b)
PoissonSolution(b::AbstractBSplineBasis, coeffs::AbstractVector) = Spline(b, coeffs)

basis(p::PoissonSolverSpline) = basis(p.quadrature)
ndofs(b::AbstractBSplineBasis) = nbasis(b)
Base.length(p::PoissonSolverSpline) = nbasis(p.quadrature)

"""
    loadvector(p::PoissonSolverSpline, f)

The Galerkin load vector ``\\int f \\phi_i \\, dx`` of the function `f`.

This is not `SimpleSplines.l2_projection`, which solves with the mass matrix as well. A Poisson
right-hand side is the load vector itself.

It returns a fresh vector, so [`solve!`](@ref) on a function allocates one per call. Build the
load vector once and solve with that where the allocation matters.
"""
function loadvector(p::PoissonSolverSpline, f)
    q = p.quadrature
    basis_values(q, 0) * (quadrature_weights(q) .* f.(quadrature_nodes(q)))
end

function solve!(result::AbstractVector, p::PoissonSolverSpline, rhs::AbstractVector)
    meanfree!(result, rhs, basis(p))
    mass_solve!(result, p.stiffness, result)
    return result
end

function solve!(result::AbstractVector, p::PoissonSolverSpline, rhs::Base.Callable)
    solve!(result, p, loadvector(p, rhs))
end

solve(p::PoissonSolverSpline, rhs::AbstractVector) = solve!(similar(rhs), p, rhs)
solve(p::PoissonSolverSpline, rhs::Base.Callable) = solve(p, loadvector(p, rhs))
