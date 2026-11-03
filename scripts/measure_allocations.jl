# Allocation and inference of the repeatedly-called paths.
#
# Run in a fresh process — a figure taken inside a session that already exercised these paths
# measures what the previous call compiled:
#
#     julia --startup-file=no --check-bounds=auto --project=. scripts/measure_allocations.jl
#
# The figures are the same under `--check-bounds=yes`, which CI forces, on Julia 1.13.3, and the
# suite's own allocation assertions are skipped only under coverage instrumentation.
#
# The claim this establishes: every solve is allocation-free, and so is evaluating a grid
# solution. Evaluating a *spline* is not: the allocation is the one `local_width` buffer
# `SimpleSplines.evaluate` takes per scalar call, reported here alongside so that the two cannot be
# confused. `SimpleSplines.evaluate_all!`, and `evaluate!` from SimpleSplines 0.4, take a
# caller-supplied buffer and allocate nothing; `Potential` holds no such buffer, by the author's
# decision (K1 in `KNOWN_ISSUES.md`).

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
    grid2 = FiniteDifferenceBasis((0.0, 1.0), NGRID; order = 2)
    grid4 = FiniteDifferenceBasis((0.0, 1.0), NGRID; order = 4)

    source(x) = 4π^2 * sin(2π * x)

    rows = Tuple{String, Int}[]
    for (name, b) in ("periodic spline" => periodic, "Dirichlet spline" => dirichlet,
        "FFT" => grid, "matrix-free order 2" => grid2, "matrix-free order 4" => grid4)
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
