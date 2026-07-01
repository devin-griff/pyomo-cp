"""CP-SAT backend.

Translates an integer-domain Pyomo model (variables, linear constraints,
``pyomo.gdp`` disjunctions, logical constraints, and a linear objective) into an
OR-Tools CP-SAT model and solves it, mapping the solution back onto the Pyomo
variables.

If the model still contains continuous variables, the backend raises and points
the user to ``TransformationFactory('cp.discretize')`` rather than discretizing
silently.

Registered as ``SolverFactory('cpsat')``. ``ortools`` is imported lazily so the
package imports without it; install ``pyomo-cp[cpsat]`` to use this backend.

Status: Phases 1-2 stub (see ROADMAP.md).
"""
from pyomo.opt import SolverFactory
from pyomo.opt.base.solvers import OptSolver


@SolverFactory.register(
    "cpsat", doc="CP-SAT (OR-Tools) backend for Pyomo (pyomo-cp)."
)
class CPSATSolver(OptSolver):
    """Solve a discretized, integer-domain Pyomo model with OR-Tools CP-SAT."""

    def __init__(self, **kwds):
        kwds.setdefault("type", "cpsat")
        super().__init__(**kwds)

    def available(self, exception_flag=False):
        try:
            import ortools  # noqa: F401
        except ImportError:
            if exception_flag:
                raise RuntimeError(
                    "The 'cpsat' backend requires OR-Tools: "
                    "pip install 'pyomo-cp[cpsat]'"
                )
            return False
        return True

    def solve(self, model, **kwds):
        raise NotImplementedError(
            "CP-SAT translation is not implemented yet (Phases 1-2; see "
            "ROADMAP.md)."
        )
