### 1D Laplace operator with periodic bc
function _apply_Δₓ!(y::AbstractVector, x::AbstractVector, h₁)
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())

    @inbounds for i in 1:nx
        i₋ = mod1(i-1, nx)
        i₊ = mod1(i+1, nx)
        y[i] = 1 / (h₁^2) * (x[i₋] - 2*x[i] + x[i₊])
    end
end

function _apply_Δₓ₄!(y::AbstractVector, x::AbstractVector, h₁) # 1D Laplace operator (4th order) with periodic bc
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())
    @inbounds for i in 1:nx
        i₋ = mod1(i-1, nx)
        i₊ = mod1(i+1, nx)
        i₋₋ = mod1(i-2, nx)
        i₊₊ = mod1(i+2, nx)
        y[i] = - 1 / (12 * h₁^2) * (x[i₊₊] - 16*x[i₊] + 30*x[i] - 16*x[i₋] + x[i₋₋])
    end
end

# L = -Δ  + R    
function _apply_Lₓ₄!(y::AbstractVector, x::AbstractVector, h₁) # 1D Laplace operator (4th order) with periodic bc
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())
    Σx = sum(x) / nx
    @inbounds for i in 1:nx
        i₋ = mod1(i-1, nx)
        i₊ = mod1(i+1, nx)
        i₋₋ = mod1(i-2, nx)
        i₊₊ = mod1(i+2, nx)
        y[i] = 1 / (12 * h₁^2) * (x[i₊₊] - 16*x[i₊] + 30*x[i] - 16*x[i₋] + x[i₋₋])
        y[i] += Σx
    end
end

### Constant Nullspace Projection
# if 1 ∈ Δ, then Δϕ = ρ is not well posed but (Δ + R)ϕ = (1 - R)ρ is. 
function _apply_Rₓ!(y::AbstractVector, x::AbstractVector) # Nullspace projection
    nx = length(x)
    length(x) == length(y) || throw(DimensionMismatch())
    Σx = sum(x) / nx
    @inbounds for i in 1:nx
        y[i] = Σx
    end
end
