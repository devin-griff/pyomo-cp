"""Discretization transformation.

Maps bounded continuous variables onto an integer grid so the model can be
solved by a finite-domain CP solver, while leaving ``pyomo.gdp`` disjunctions
intact (only the variables are changed, not the disjunctive structure).

This is a general, solver-agnostic transformation and an *explicit* step: CP
backends do not run it automatically, because discretizing is a modelling
decision (the integrality assumption) that changes the problem.

Registered as ``TransformationFactory('cp.discretize')``.

Two modes:

* ``step == 1`` (unit grid): each bounded continuous variable is set to the
  integer domain in place, with bounds tightened inward to the enclosed integer
  range. For integer-valued bounds this is exact; the solution lands directly in
  the original variables.
* ``step != 1``: each bounded continuous variable ``x in [lb, ub]`` is replaced
  by ``x = lb + step * x_int`` with an integer grid variable
  ``x_int in [0, floor((ub-lb)/step)]``, substituted throughout the constraints
  and objective. The original variable is fixed and its value is recovered from
  the grid variable after the solve (the ``cpsat`` backend does this
  automatically via the map stored on the model).
"""
import math

from pyomo.core import (
    Block,
    Constraint,
    Integers,
    Objective,
    Transformation,
    TransformationFactory,
    Var,
    value,
)
from pyomo.core.expr import replace_expressions
from pyomo.repn import generate_standard_repn

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
        Grid resolution (default 1, a unit grid).
    """

    def _apply_to(self, model, step=1, **kwds):
        if step <= 0:
            raise ValueError("pyomo-cp: cp.discretize step must be positive.")
        if step == 1:
            self._unit_grid(model)
        else:
            self._general_grid(model, float(step))

    @staticmethod
    def _single_var_pin(c):
        """If equality constraint ``c`` reduces to ``coef*x + const == rhs`` in a
        single variable ``x``, return ``(x, implied_value)``; else ``None``."""
        if not c.equality:
            return None
        repn = generate_standard_repn(c.body)
        if not repn.is_linear() or len(repn.linear_vars) != 1:
            return None
        coef = float(repn.linear_coefs[0])
        if coef == 0:
            return None
        rhs = float(value(c.upper))
        return repn.linear_vars[0], (rhs - float(repn.constant)) / coef

    def _check_offgrid_pins(self, model, lb_of, step, n_of):
        """Raise if any single-variable equality pins a discretized variable to a
        value off its grid. This is the common, cheaply detectable way
        discretization turns a feasible model infeasible (e.g. a dimension fixed
        to an odd value on an even-step grid), caught here rather than surfacing
        as a puzzling ``infeasible`` at solve time. General (multi-variable)
        infeasibility is not detectable without solving and is not checked."""
        for c in model.component_data_objects(
            Constraint, active=True, descend_into=True
        ):
            pin = self._single_var_pin(c)
            if pin is None:
                continue
            v, val = pin
            if id(v) not in lb_of:
                continue
            lb, n = lb_of[id(v)], n_of[id(v)]
            k = (val - lb) / step
            kr = round(k)
            if abs(k - kr) > _TOL or kr < 0 or kr > n:
                gmax = lb + step * n
                raise ValueError(
                    f"pyomo-cp: constraint '{c.name}' pins variable '{v.name}' to "
                    f"{val:g}, which is not on its discretization grid "
                    f"({lb:g}, {lb + step:g}, ..., {gmax:g} at step {step:g}). "
                    f"Discretizing would make the model infeasible. Use a step "
                    f"that divides the required values (a unit or fractional "
                    f"grid), or adjust the data/bounds."
                )

    @staticmethod
    def _bounds_or_raise(v):
        lb, ub = v.bounds
        if lb is None or ub is None:
            raise ValueError(
                f"pyomo-cp: cannot discretize unbounded variable '{v.name}'; "
                f"give it finite bounds first."
            )
        return lb, ub

    def _unit_grid(self, model):
        lb_of, n_of = {}, {}
        for v in model.component_data_objects(Var, active=True, descend_into=True):
            if v.fixed or not v.is_continuous():
                continue
            lb, ub = self._bounds_or_raise(v)
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
            lb_of[id(v)], n_of[id(v)] = ilb, iub - ilb
        self._check_offgrid_pins(model, lb_of, 1, n_of)

    def _general_grid(self, model, step):
        cont = [
            v
            for v in model.component_data_objects(Var, active=True, descend_into=True)
            if v.is_continuous() and not v.fixed
        ]
        if not cont:
            return
        if not hasattr(model, "_cp_disc"):
            model._cp_disc = Block()

        sub = {}
        disc = []
        lb_of, n_of = {}, {}
        for x in cont:
            lb, ub = self._bounds_or_raise(x)
            n = math.floor((ub - lb) / step + _TOL)
            if n < 0:
                raise ValueError(
                    f"pyomo-cp: variable '{x.name}' has no grid point at step "
                    f"{step} in [{lb}, {ub}]."
                )
            xi = Var(domain=Integers, bounds=(0, n))
            model._cp_disc.add_component(f"v{len(disc)}", xi)
            sub[id(x)] = lb + step * xi
            disc.append((x, xi, lb, step))
            lb_of[id(x)], n_of[id(x)] = lb, n

        # Catch off-grid pins before rewriting the constraints (the check reads
        # the original single-variable equalities).
        self._check_offgrid_pins(model, lb_of, step, n_of)

        for c in model.component_data_objects(Constraint, active=True, descend_into=True):
            c.set_value(replace_expressions(c.expr, sub))
        for o in model.component_data_objects(Objective, active=True):
            o.set_value(replace_expressions(o.expr, sub))

        for x, xi, lb, stp in disc:
            x.fix(float(lb))  # placeholder; real value recovered on descale
        model._pyomo_cp_disc = getattr(model, "_pyomo_cp_disc", []) + disc
