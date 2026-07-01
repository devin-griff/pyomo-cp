"""Discretization transformation.

Maps each bounded continuous variable ``x in [lb, ub]`` onto an integer variable
``x_int`` via ``x = lb + step * x_int``, and scales the referencing constraints
and objective to integer coefficients, so the model can be solved by a
finite-domain CP solver. ``pyomo.gdp`` disjunctions are preserved (the variables
inside disjuncts are discretized, but the disjunctions are NOT reformulated to
big-M/hull), so a CP backend can map them to native reified constraints.

This is deliberately a *general, solver-agnostic* transformation and an
*explicit* step. CP backends do not run it automatically: discretizing is a
modelling decision (the integrality assumption) that changes the problem, so the
user should choose the resolution and apply it knowingly.

Registered as ``TransformationFactory('cp.discretize')``.

Status: Phase 3 stub (see ROADMAP.md).
"""
from pyomo.core import Transformation, TransformationFactory


@TransformationFactory.register(
    "cp.discretize",
    doc="Discretize bounded continuous variables onto integer grids (pyomo-cp).",
)
class DiscretizeTransformation(Transformation):
    """Discretize bounded continuous variables onto integer grids.

    Parameters (planned)
    ---------------------
    step : float or mapping
        Grid resolution; a global default or a per-variable mapping.

    Exactness: for integer-data models the transformation is exact; for models
    with non-integer coefficients it is an approximation controlled by ``step``.
    """

    def _apply_to(self, model, step=1.0, **kwds):
        raise NotImplementedError(
            "cp.discretize is not implemented yet (Phase 3; see ROADMAP.md)."
        )
