using PoissonSolvers
using PoissonSolvers: discretise, ndofs, rhs, update!
using SimpleSplines
using Test

# One `Potential` per backend, exercised identically. Everything asserted here is backend
# independent, so anything that holds for one and not another is a defect in that one.
#
# The evaluation tolerance is the exception, and is carried per backend for a reason rather than
# loosened to whichever is worse. A spline evaluates where it is asked; the grid solution snaps to
# the nearest point, so its error is O(Δx |φ'|) ≈ 0.05 here whatever the solve does. Giving the
# spline the grid's tolerance would hide a real regression in it.
const BASES = ("spline" => (PeriodicBasisSpline((0.0, 1.0), 5, 32), 1e-3),
    "grid" => (FFTWBasis((0.0, 1.0), 64), 0.06),
    "matrix-free" => (FiniteDifferenceBasis((0.0, 1.0), 64; order = 4), 0.06))

source(x) = 4π^2 * sin(2π * x)
exact(x) = sin(2π * x)

@testset "$name" for (name, (b, atol)) in BASES
    potential = Potential(b, source)

    @testset "construction" begin
        @test basis(potential) === b
        @test length(coefficients(potential)) == ndofs(b)
        @test length(rhs(potential)) == ndofs(b)
        @test all(iszero, rhs(potential))

        # The single-argument form solves with a zero source, so the potential is zero.
        @test Potential(b)(0.3) == 0.0
    end

    @testset "evaluation" begin
        xs = 0.05:0.05:0.95
        @test maximum(abs, potential.(xs) .- exact.(xs)) < atol
        @test potential(0.3) == potential.potential(0.3)
        @test potential(collect(xs)) == potential.(xs)
        @test potential(collect(xs), 1) == potential.(xs, 1)
        @test derivative(potential).(xs) == potential.(xs, 1)
    end

    @testset "update! writes through the solution" begin
        # The solution shares its coefficient array rather than copying it, so an in-place
        # re-solve is visible through the functor with nothing rebuilt. This is the contract the
        # whole `Potential` type exists for, and it is invisible from the outside.
        coeffs = coefficients(potential)
        @test update!(potential, zero(coeffs)) === potential
        @test coefficients(potential) === coeffs
        @test all(iszero, coeffs)
        @test potential(0.3) == 0.0

        # A function is accepted too, and each backend reduces it to its own right-hand side.
        # Restoring the source restores the solution, so `update!` is not one-way.
        update!(potential, source)
        @test potential(0.3) ≈ exact(0.3) atol=atol

        # With no argument it re-solves from the stored right-hand side, which is zero.
        update!(potential)
        @test all(iszero, coefficients(potential))
    end
end

@testset "the element type follows the basis" begin
    # The zero right-hand side of the single-argument constructor was Float64 whatever the basis
    # was, and the two backends failed differently. The grid basis threw a MethodError from `mul!`,
    # because the transform is built for the basis element type. The spline basis took the wider
    # vector and returned a Float64 potential without complaint, which is why both are asserted
    # here: the loud failure is the easy one to catch.
    for b in (FFTWBasis(Float32.((0.0, 1.0)), 64),
        FiniteDifferenceBasis(Float32.((0.0, 1.0)), 64),
        DirichletBasisSpline(Float32.((0.0, 1.0)), 5, 32))
        @test eltype(b) == Float32

        potential = Potential(b)
        @test eltype(rhs(potential)) == Float32
        @test eltype(coefficients(potential)) == Float32
        @test potential(0.3f0) == 0
    end
end

# The reading is taken through a barrier whose arguments have concrete types, after a warm call,
# so that it measures the call itself rather than its compilation or a dynamic dispatch around it.
allocations(f::F, a::A) where {F, A} = (f(a); @allocated f(a))
allocations(f::F, a::A, b::B) where {F, A, B} = (f(a, b); @allocated f(a, b))

# `derivative(p, 1)(x)` as one call, so that `@inferred` sees the type that `derivative` returns.
derivative_at(p, x) = derivative(p, 1)(x)

# Skipped where a run forces --check-bounds=yes or instruments coverage.
measured() = Base.JLOptions().check_bounds == 0 && Base.JLOptions().code_coverage == 0

@testset "update! is inferred and allocates nothing, $T" for T in (Float64, Float32)
    domain = T.((0.0, 1.0))
    for b in (FFTWBasis(domain, 64), FiniteDifferenceBasis(domain, 64; order = 4),
        PeriodicBasisSpline(domain, 5, 32), DirichletBasisSpline(domain, 5, 32))
        potential = Potential(b)
        ρ = discretise(potential.solver, source)
        @test @inferred(update!(potential, ρ)) === potential
        @test !all(iszero, coefficients(potential))

        if measured()
            @test allocations(update!, potential, ρ) == 0
        end
        @test allocations((p, ρ) -> similar(ρ), potential, ρ) > 0
    end
end

@testset "a grid potential evaluates without allocating, $T" for T in (Float64, Float32)
    domain = T.((0.0, 1.0))
    x = T(0.3)
    fft = Potential(FFTWBasis(domain, 64), source)
    fd = Potential(FiniteDifferenceBasis(domain, 64; order = 4), source)

    @test @inferred(fft(x, 1)) isa T
    @test @inferred(derivative_at(fft, x)) isa T
    @test @inferred(fd(x)) isa T

    if measured()
        @test allocations(fft, x, 1) == 0
        @test allocations(derivative_at, fft, x) == 0
        @test allocations(fd, x) == 0
    end
    @test allocations((p, x) -> similar(coefficients(p)), fft, x) > 0
    @test allocations(x -> similar(coefficients(fd)), x) > 0
end

@testset "a spline potential evaluates with an inferred type, $T" for T in (Float64, Float32)
    # A spline solution is a SimpleSplines `Spline`, whose point evaluation allocates a buffer per
    # call, so only the inferred half is asserted here.
    domain = T.((0.0, 1.0))
    x = T(0.3)
    for b in (PeriodicBasisSpline(domain, 5, 32), DirichletBasisSpline(domain, 5, 32))
        spline = Potential(b, source)
        @test @inferred(spline(x)) isa T
        @test @inferred(spline(x, 1)) isa T
        @test @inferred(derivative_at(spline, x)) isa T
    end
end
