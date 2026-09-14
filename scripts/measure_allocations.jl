# Allocation and inference of the repeatedly-called paths.
#
# Run in a fresh process — a figure taken inside a session that already exercised these paths
# measures what the previous call compiled:
#
#     julia --startup-file=no --check-bounds=auto --project=. scripts/measure_allocations.jl
#
# `--check-bounds=auto` matters. `Pkg.test()` forces `--check-bounds=yes` up to Julia 1.12, which
# inflates every figure here, and is why the suite's own allocation assertions are guarded on
# `Base.JLOptions().check_bounds == 0` rather than made unconditionally.
#
# The claim this establishes: every solve is allocation-free, and so is evaluating a grid
# solution. Evaluating a *spline* is not, and the residue is not this package's — it is the one
# `local_width` buffer `SimpleSplines.evaluate` takes per scalar call, reported here alongside so
# that the two cannot be confused. `SimpleSplines.evaluate_all!`, which takes a caller-supplied
# buffer, allocates nothing, so the fix is upstream rather than a scratch buffer here.

using PoissonSolvers
using PoissonSolvers: update!
using Printf
using SimpleSplines
using Test

const ORDER = 5
const NCELLS = 32
const NGRID = 64

"""Run `f` a few times, then return its allocation in bytes."""
function allocation(f)
    for _ in 1:3
        f()
    end
    @allocated f()
end

function report()
    periodic = PeriodicBasisSpline((0.0, 1.0), ORDER, NCELLS)
    dirichlet = DirichletBasisSpline((0.0, 1.0), ORDER, NCELLS)
    grid = FFTWBasis((0.0, 1.0), NGRID)

    source(x) = 4π^2 * sin(2π * x)

    rows = Tuple{String, Int}[]
    for (name, b) in ("periodic spline" => periodic, "Dirichlet spline" => dirichlet,
        "grid" => grid)
        solver = PoissonSolver(b)
        ρ = rand(length(solver))
        φ = similar(ρ)
        potential = Potential(b, source)

        push!(rows, ("solve!  $(name)", allocation(() -> solve!(φ, solver, ρ))))
        push!(rows, ("update! $(name)", allocation(() -> update!(potential, ρ))))
        push!(rows, ("p(x)    $(name)", allocation(() -> potential(0.3))))
        push!(rows, ("p(x, 1) $(name)", allocation(() -> potential(0.3, 1))))
    end

    û = rand(nbasis(periodic))
    buffer = zeros(SimpleSplines.local_width(periodic))
    push!(rows, (
        "SimpleSplines.evaluate(b, û, x)", allocation(() -> evaluate(periodic, û, 0.3))))
    push!(rows,
        ("SimpleSplines.evaluate_all!(buf, b, x)",
            allocation(() -> evaluate_all!(buffer, periodic, 0.3, 0))))

    println("check_bounds = ", Base.JLOptions().check_bounds, "  (0 = auto, 1 = yes, 2 = no)\n")
    for (name, bytes) in rows
        @printf("%-40s %6d bytes\n", name, bytes)
    end
end

report()
