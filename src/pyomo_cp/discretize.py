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
from fractions import Fraction

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
    def _resolve_grid(v, lb_of, step):
        """Return ``(lb, step)`` describing the grid ``v = lb + step*k`` (integer
        k) that ``v`` will live on, or ``None`` if it can't be determined (an
        already-integer variable contributes a unit grid at its lower bound; a
        continuous variable that isn't being discretized can't be reasoned about,
        so the constraint is skipped)."""
        if id(v) in lb_of:
            return lb_of[id(v)], step
        if v.is_integer() or v.is_binary():
            lb = v.lb
            if lb is not None:
                return float(lb), 1.0
        return None

    def _check_offgrid_pins(self, model, lb_of, step, n_of):
        """Raise if a linear equality has no solution on the discretization grid.

        Substituting ``x = lb + step*k`` (integer k) turns each equality
        ``sum a_i x_i == rhs`` into ``sum (a_i step_i) k_i == R``. Cleared to
        integer coefficients, that has an integer solution only if ``gcd`` of the
        coefficients divides the right-hand side. When it does not, no grid
        assignment can satisfy the constraint, so discretizing would make the
        model infeasible -- caught here rather than surfacing as a puzzling
        ``infeasible`` at solve time.

        The test ignores variable bounds, so it is a *necessary* condition: it
        only ever rejects genuinely-infeasible discretizations, never good ones,
        and it stays silent when a solution might exist (bound-tightness and
        multi-constraint interactions are left to the solver)."""
        for c in model.component_data_objects(
            Constraint, active=True, descend_into=True
        ):
            if not c.equality:
                continue
            repn = generate_standard_repn(c.body)
            if not repn.is_linear() or not repn.linear_vars:
                continue
            coefs = [float(a) for a in repn.linear_coefs]
            grids = [self._resolve_grid(v, lb_of, step) for v in repn.linear_vars]
            if any(g is None for g in grids):
                continue  # a variable we can't reason about; defer to the solver

            rhs = float(value(c.upper)) - float(repn.constant)
            # sum a_i (lb_i + step_i k_i) == rhs  =>  sum (a_i step_i) k_i == R
            residual = rhs - sum(a * lb for a, (lb, st) in zip(coefs, grids))
            terms = [a * st for a, (lb, st) in zip(coefs, grids)]

            fracs = [Fraction(x).limit_denominator(10 ** 6) for x in terms + [residual]]
            denom = 1
            for f in fracs:
                denom = denom * f.denominator // math.gcd(denom, f.denominator)
            ints = [int(f * denom) for f in fracs]
            coeff_gcd = 0
            for b in ints[:-1]:
                coeff_gcd = math.gcd(coeff_gcd, abs(b))
            target = ints[-1]

            feasible = (target == 0) if coeff_gcd == 0 else (target % coeff_gcd == 0)
            if not feasible:
                raise ValueError(self._offgrid_message(c, repn, coefs, grids, n_of))

    @staticmethod
    def _offgrid_message(c, repn, coefs, grids, n_of):
        if len(repn.linear_vars) == 1:
            v = repn.linear_vars[0]
            lb, st = grids[0]
            val = (float(value(c.upper)) - float(repn.constant)) / coefs[0]
            n = n_of.get(id(v))
            if n is not None:
                grid = f"({lb:g}, {lb + st:g}, ..., {lb + st * n:g} at step {st:g})"
            else:
                grid = f"(step {st:g} from {lb:g})"
            return (
                f"pyomo-cp: constraint '{c.name}' pins variable '{v.name}' to "
                f"{val:g}, which is not on its discretization grid {grid}. "
                f"Discretizing would make the model infeasible. Use a step that "
                f"divides the required values (a unit or fractional grid), or "
                f"adjust the data/bounds."
            )
        return (
            f"pyomo-cp: equality constraint '{c.name}' has no solution on the "
            f"discretization grid (no integer grid assignment satisfies it), so "
            f"discretizing would make the model infeasible. Use a step that "
            f"divides the constraint's data, or adjust it."
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
