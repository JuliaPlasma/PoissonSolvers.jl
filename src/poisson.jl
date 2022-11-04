
abstract type PoissonSolver{dType} end


"""
Solves the Poisson equation for periodic boundary conditions

```
function solve!(p::PoissonSolver, x::AbstractVector, w::AbstractVector) end
```

"""
function solve! end


"""
Evaluates the charge density field
"""
function eval_density end

"""
Evaluates the electrostatic potential
"""
function eval_potential end

"""
Evaluates the electric field
"""
function eval_field end



function eval_density(p::PoissonSolver{DT}, x::AbstractVector{DT}) where {DT}
    [eval_density(p,x_) for x_ in x]
end

function eval_potential(p::PoissonSolver{DT}, x::AbstractVector{DT}) where {DT}
    [eval_potential(p,x_) for x_ in x]
end

function eval_field(p::PoissonSolver{DT}, x::AbstractVector{DT}) where {DT}
    [eval_field(p,x_) for x_ in x]
end


function eval_density!(y::AbstractVector{DT}, p::PoissonSolver{DT}, x::AbstractVector{DT}) where {DT}
    for i in eachindex(x,y)
        y[i] = eval_density(p,x[i])
    end
    return y
end

function eval_potential!(y::AbstractVector{DT}, p::PoissonSolver{DT}, x::AbstractVector{DT}) where {DT}
    for i in eachindex(x,y)
        y[i] = eval_potential(p,x[i])
    end
    return y
end

function eval_field!(y::AbstractVector, p::PoissonSolver, x::AbstractVector)
    for i in eachindex(x,y)
        y[i] = eval_field(p,x[i])
    end
    return y
end

function eval_field!(y::AbstractMatrix, p::PoissonSolver, x::AbstractMatrix)
    @assert size(x, 1) == 1
    @assert size(y, 1) == 1

    X = reshape(x, size(x,2))
    Y = reshape(y, size(y,2))

    eval_field!(Y, p, X)
end
