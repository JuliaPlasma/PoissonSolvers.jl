
abstract type PoissonSolver{dType} end


"""
Solves the Poisson equation for periodic boundary conditions

```
function solve!(p::PoissonSolver, x::AbstractVector, w::AbstractVector) end
```

"""
function solve! end



"""
Takes a basis and a coefficient vector and returns an appropriate object 
that allows to evaluate the corresponding solution via a functor via

```
sol = PoissonSolution(basis, coeffs)
sol(x)
```

"""
function PoissonSolution end


function evalsolution(basis, coeffs, x::AbstractVector)
    [evalsolution(basis, coeffs, x_) for x_ in x]
end
