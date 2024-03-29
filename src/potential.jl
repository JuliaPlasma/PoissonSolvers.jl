
struct Potential{BT, CT <: AbstractVector, PT, ST <: PoissonSolver}
    basis::BT
    coefficients::CT
    potential::PT
    solver::ST

    function Potential(basis, init_rhs) where {DT}
        solver = PoissonSolver(basis)
        coeffs = solve(solver, init_rhs)
        potential = PoissonSolution(basis, coeffs)

        new{typeof(basis), typeof(coeffs), typeof(potential), typeof(solver)}(basis, coeffs, potential, solver)
    end
end

Potential(basis) = Potential(basis, zeros(length(basis)))

(p::Potential)(x::Number) = p.potential(x)
(p::Potential)(x::AbstractArray) = p.potential.(x)

(p::Potential)(x::Number, derivative::Derivative) = (derivative * p.potential)(x)
(p::Potential)(x::AbstractArray, derivative::Derivative) = (derivative * p.potential).(x)


function update!(p::Potential, rhs)
    solve!(p.coefficients, p.solver, rhs)
    return p
end
