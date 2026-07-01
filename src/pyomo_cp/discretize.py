"""Discretization transformation.

Maps bounded continuous variables onto an integer grid so the model can be
solved by a finite-domain CP solver, while leaving ``pyomo.gdp`` disjunctions
intact (only the variables are changed, not the disjunctive structure).

This is a general, solver-agnostic transformation and an *explicit* step: CP
backends do not run it automatically, because discretizing is a modelling
decision (the integrality assumption) that changes the problem, so the user
should choose the resolution and apply it knowingly.

Registered as ``TransformationFactory('cp.discretize')``.

Current scope: unit grid (``step=1``). Each bounded continuous variable is set
to the integer domain, with its bounds tightened inward to the enclosed integer
range. For integer-valued bounds this is exact; the solution lands directly in
the original variables (no descaling needed). General non-unit grids are a
planned extension (see ROADMAP.md).
"""
import math

from pyomo.core import Integers, Var, Transformation, TransformationFactory

_TOL = 1e-9


@TransformationFactory.register(
    "cp.discretize",
    doc="Discretize bounded continuous variables onto integer grids (pyomo-cp).",
)
class DiscretizeTransformation(Transformation):
    """Discretize bounded continuous variables onto an integer grid.

    Parameters
    ----------
    step : float
        Grid resolution. Only ``step == 1`` (unit grid) is supported so far.
    """

    def _apply_to(self, model, step=1, **kwds):
        if step != 1:
            raise NotImplementedError(
                "pyomo-cp: cp.discretize currently supports only step=1 (unit "
                "grid); general step sizes are planned (see ROADMAP.md)."
            )
        for v in model.component_data_objects(Var, active=True, descend_into=True):
            if v.fixed or not v.is_continuous():
                continue
            lb, ub = v.bounds
            if lb is None or ub is None:
                raise ValueError(
                    f"pyomo-cp: cannot discretize unbounded variable '{v.name}'; "
                    f"give it finite bounds first."
                )
            ilb = math.ceil(lb - _TOL)
            iub = math.floor(ub + _TOL)
            if ilb > iub:
                raise ValueError(
                    f"pyomo-cp: variable '{v.name}' has no integer point in its "
                    f"bounds [{lb}, {ub}]."
                )
            v.domain = Integers
            v.setlb(ilb)
            v.setub(iub)
