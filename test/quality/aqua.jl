using Aqua
using PoissonSolvers
using Test

# Package-level quality assurance: type piracy, method ambiguities, stale and duplicated
# dependencies, undefined exports, unbound type parameters, `Project.toml` validity. These are
# the faults the rest of the suite is structurally unable to see — it exercises behaviour, and
# every one of these is a property of the package as a whole.
#
# Two of them are live concerns here. `undefined_exports` guards the export list, which names
# backend entry points that are defined in different files from the module that exports them, so
# a rename in one place and not the other is invisible until a caller reaches for the name.
# `stale_deps` guards the dependency list against exactly the residue a backend swap leaves
# behind.
#
# Piracy is the third. The solvers extend `SimpleSplines`' `basis`, `coefficients` and
# `derivative`, and `Base.length`, on types this package owns — which is extension, not piracy.
# A method that owned neither the function nor an argument type would be, and this is what says
# so.
Aqua.test_all(PoissonSolvers)
