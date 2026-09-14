@doc raw"""
    Potential(basis)
    Potential(basis, rhs)

The solution of ``-\Delta \phi = \rho`` on `basis`, kept together with the solver that produced it
so that it can be re-solved in place as `rhs` changes.

Callable: `p(x)` evaluates ``\phi``, and `p(x, d)` its `d`-th derivative. `basis(p)` and
`coefficients(p)` reach the basis and the coefficient vector; both live inside the solution rather
than being stored a second time.

`update!(p, rhs)` overwrites the coefficients in place. That is visible through `p` itself because
the solution **shares** its coefficient array rather than copying it — the contract every backend
here keeps, and what makes a re-solve from a vector free of allocation.
"""
struct Potential{PT, ST <: PoissonSolver, CT <: AbstractVector}
    potential::PT
    solver::ST
    rhs::CT

    function Potential(b, init_rhs)
        solver = PoissonSolver(b)
        coeffs = solve(solver, init_rhs)
        potential = PoissonSolution(b, coeffs)
        new{typeof(potential), typeof(solver), typeof(coeffs)}(potential, solver, zero(coeffs))
    end
end

Potential(b) = Potential(b, zeros(ndofs(b)))

basis(p::Potential) = basis(p.potential)
coefficients(p::Potential) = coefficients(p.potential)

(p::Potential)(x::Number) = p.potential(x)
(p::Potential)(x::Number, d::Integer) = p.potential(x, d)
(p::Potential)(x::AbstractArray) = p.potential.(x)
(p::Potential)(x::AbstractArray, d::Integer) = derivative(p, d).(x)

derivative(p::Potential, d = 1) = derivative(p.potential, d)

"""
    rhs(p::Potential)

The stored right-hand side buffer of `p`, which [`update!`](@ref) solves with by default.

It is a scratch array for the caller to fill, and starts out zero: the right-hand side given to
the constructor is solved with, not stored here. So `update!(p)` on a freshly built `p` solves
for a zero source.
"""
rhs(p::Potential) = p.rhs

"""
    update!(p::Potential, rhs = rhs(p))

Re-solve `p` for the right-hand side `rhs`, in place, and return `p`.

`rhs` is either a vector or a function, as for [`solve!`](@ref). The new coefficients are written
into the array the solution already holds, so the result is visible through `p` without anything
being rebuilt.

With a vector — the stored [`rhs`](@ref) buffer in particular — nothing is allocated. A function
is sampled onto a fresh vector first, which allocates one per call.
"""
function update!(p::Potential, rhs = p.rhs)
    solve!(coefficients(p), p.solver, rhs)
    return p
end
