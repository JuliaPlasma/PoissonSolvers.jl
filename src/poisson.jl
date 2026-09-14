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
