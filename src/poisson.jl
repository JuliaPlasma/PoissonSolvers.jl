abstract type PoissonSolver{dType} end

"""
    PoissonSolver(basis)

The solver this `basis` selects for ``-\\Delta \\phi = \\rho``.
"""
function PoissonSolver end

"""
    PoissonSolution(basis, coefficients)

A callable solution built from `basis` and `coefficients`, so that `sol(x)` evaluates it.

The result **shares** `coefficients` rather than copying it, which is what lets
[`update!`](@ref) re-solve in place.
"""
function PoissonSolution end

"""
    ndofs(basis)

The number of degrees of freedom of `basis`.

The backends spell this differently — `nbasis` for a spline basis, `length` for a grid — and
[`Potential`](@ref) should not have to know which one it holds.
"""
function ndofs end

"""
    solve!(coefficients, solver, rhs)
    solve(solver, rhs)

Solve ``-\\Delta \\phi = \\rho`` for the coefficients of ``\\phi``.

`rhs` is either the discrete right-hand side or a function, which each backend reduces to one.
"""
function solve! end

@doc (@doc solve!) function solve end

"""
    discretise(solver, f)

The discrete right-hand side of the function `f` on `solver`'s basis: for a grid, the samples the
solver takes; for a spline, the Galerkin load vector.
"""
function discretise end

# A wrong length otherwise reaches a planned transform or a factorisation, which reports it
# against a plan the caller never made. Naming the three lengths says what disagrees.
function checklengths(result, rhs, p::PoissonSolver)
    length(result) == length(rhs) == length(p) || throw(DimensionMismatch(
        "the solver has $(length(p)) degrees of freedom, but the right-hand side has " *
        "$(length(rhs)) and the result $(length(result))"))
    return nothing
end

function solve!(coeffs::AbstractVector, p::PoissonSolver, rhs::Base.Callable)
    solve!(coeffs, p, discretise(p, rhs))
end
solve(p::PoissonSolver, rhs::AbstractVector) = solve!(similar(rhs), p, rhs)
solve(p::PoissonSolver, rhs::Base.Callable) = solve(p, discretise(p, rhs))
