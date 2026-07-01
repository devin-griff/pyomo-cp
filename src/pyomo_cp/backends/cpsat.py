"""CP-SAT backend.

Translates a Pyomo model (variables, linear constraints, ``pyomo.gdp``
disjunctions, logical constraints, and a linear objective) into an OR-Tools
CP-SAT model, solves it, and loads the solution back onto the Pyomo variables.

Disjunctions map to native reified constraints; logical constraints and Boolean
variables map to CP-SAT's boolean constraints. Continuous variables are rejected
with a pointer to ``TransformationFactory('cp.discretize')``.

Registered as ``SolverFactory('cpsat')``. ``ortools`` is imported lazily so the
package imports without it; install ``pyomo-cp[cpsat]`` to use this backend.
"""
from fractions import Fraction
from math import gcd

from pyomo.core import (
    BooleanVar,
    Block,
    Constraint,
    Objective,
    Var,
    minimize,
    value,
)
from pyomo.core.base.boolean_var import BooleanVarData
from pyomo.gdp import Disjunction
from pyomo.opt import SolverFactory, SolverResults, SolverStatus, TerminationCondition
from pyomo.opt.base.solvers import OptSolver
from pyomo.repn import generate_standard_repn

try:
    from pyomo.core import LogicalConstraint
except ImportError:  # pragma: no cover
    LogicalConstraint = None

_INT_TOL = 1e-9


def _as_int(x, what="value"):
    xi = round(x)
    if abs(x - xi) > _INT_TOL:
        raise ValueError(
            f"pyomo-cp: the CP-SAT backend needs an integer {what}, got {x}. "
            f"Scale the data or apply TransformationFactory('cp.discretize')."
        )
    return int(xi)


def _int_scale(values, max_denom=10**6):
    """Smallest positive integer s such that s*v is integral for every v."""
    lcm = 1
    for v in values:
        d = Fraction(v).limit_denominator(max_denom).denominator
        lcm = lcm * d // gcd(lcm, d)
    return lcm


def _linear_expr(repn, varmap):
    """CP-SAT linear expression from a Pyomo standard repn (integer coefficients
    required; used for the objective)."""
    expr = _as_int(repn.constant, "constant")
    for coef, var in zip(repn.linear_coefs, repn.linear_vars):
        if id(var) in varmap:
            expr = expr + _as_int(coef, "coefficient") * varmap[id(var)][1]
        else:
            expr = expr + _as_int(coef * value(var), "coefficient")
    return expr


def _emit_constraint(cpm, c, varmap, enforce):
    """Add a linear constraint, scaling coefficients to integers and reifying on
    the enforcing indicators if any."""
    repn = generate_standard_repn(c.body)
    if not repn.is_linear():
        raise ValueError(
            f"pyomo-cp: constraint '{c.name}' is nonlinear; the CP-SAT backend "
            f"supports linear constraints only."
        )
    coefs = [float(x) for x in repn.linear_coefs]
    const = float(repn.constant)
    lower = None if c.lower is None else float(value(c.lower))
    upper = None if c.upper is None else float(value(c.upper))
    s = _int_scale(coefs + [const] + [x for x in (lower, upper) if x is not None])

    expr = int(round(const * s))
    for coef, var in zip(repn.linear_coefs, repn.linear_vars):
        if id(var) in varmap:
            expr = expr + int(round(coef * s)) * varmap[id(var)][1]
        else:
            expr = expr + int(round(coef * value(var) * s))

    lits = list(enforce)

    def add(relation):
        ct = cpm.Add(relation)
        if lits:
            ct.OnlyEnforceIf(lits)

    if c.equality:
        add(expr == int(round(lower * s)))
    else:
        if lower is not None:
            add(expr >= int(round(lower * s)))
        if upper is not None:
            add(expr <= int(round(upper * s)))


# --- boolean / logical translation ------------------------------------------

def _get_bool(cpm, bv, boolmap, varmap):
    if id(bv) in boolmap:
        return boolmap[id(bv)]
    b = cpm.NewBoolVar(bv.name)
    boolmap[id(bv)] = b
    try:  # link to an associated binary variable if one exists and is in scope
        assoc = bv.get_associated_binary()
    except Exception:  # noqa: BLE001
        assoc = None
    if assoc is not None and id(assoc) in varmap:
        cpm.Add(b == varmap[id(assoc)][1])
    return b


def _to_literal(cpm, expr, boolmap, varmap):
    """Return a CP-SAT literal equivalent to a boolean expression, adding
    reification as needed. Counting expressions (Exactly/AtMost/AtLeast) are
    supported only at the top level of a LogicalConstraint, not nested."""
    from pyomo.core.expr.logical_expr import (
        AndExpression,
        EquivalenceExpression,
        ImplicationExpression,
        NotExpression,
        OrExpression,
    )

    if isinstance(expr, BooleanVarData):
        return _get_bool(cpm, expr, boolmap, varmap)
    if expr is True or expr is False:
        b = cpm.NewBoolVar("_const")
        cpm.Add(b == (1 if expr else 0))
        return b
    if isinstance(expr, NotExpression):
        return _to_literal(cpm, expr.args[0], boolmap, varmap).Not()
    if isinstance(expr, AndExpression):
        lits = [_to_literal(cpm, a, boolmap, varmap) for a in expr.args]
        aux = cpm.NewBoolVar("_and")
        cpm.AddMinEquality(aux, lits)
        return aux
    if isinstance(expr, OrExpression):
        lits = [_to_literal(cpm, a, boolmap, varmap) for a in expr.args]
        aux = cpm.NewBoolVar("_or")
        cpm.AddMaxEquality(aux, lits)
        return aux
    if isinstance(expr, ImplicationExpression):  # a -> b  ==  (not a) or b
        a = _to_literal(cpm, expr.args[0], boolmap, varmap)
        b = _to_literal(cpm, expr.args[1], boolmap, varmap)
        aux = cpm.NewBoolVar("_impl")
        cpm.AddMaxEquality(aux, [a.Not(), b])
        return aux
    if isinstance(expr, EquivalenceExpression):  # a <-> b
        a = _to_literal(cpm, expr.args[0], boolmap, varmap)
        b = _to_literal(cpm, expr.args[1], boolmap, varmap)
        aux = cpm.NewBoolVar("_equiv")
        cpm.Add(a == b).OnlyEnforceIf(aux)
        cpm.Add(a + b == 1).OnlyEnforceIf(aux.Not())
        return aux
    raise NotImplementedError(
        f"pyomo-cp: boolean expression '{type(expr).__name__}' is not supported "
        f"as a sub-expression (counting expressions only at the top level)."
    )


def _emit_logical(cpm, lc, boolmap, varmap, enforce):
    from pyomo.core.expr.logical_expr import (
        AtLeastExpression,
        AtMostExpression,
        ExactlyExpression,
    )

    expr = lc.expr
    lits_enforce = list(enforce)
    if isinstance(expr, (AtLeastExpression, AtMostExpression, ExactlyExpression)):
        n = int(value(expr.args[0]))
        lits = [_to_literal(cpm, a, boolmap, varmap) for a in expr.args[1:]]
        total = sum(lits)
        if isinstance(expr, AtLeastExpression):
            ct = cpm.Add(total >= n)
        elif isinstance(expr, AtMostExpression):
            ct = cpm.Add(total <= n)
        else:
            ct = cpm.Add(total == n)
        if lits_enforce:
            ct.OnlyEnforceIf(lits_enforce)
    else:
        lit = _to_literal(cpm, expr, boolmap, varmap)
        ct = cpm.AddBoolAnd([lit])
        if lits_enforce:
            ct.OnlyEnforceIf(lits_enforce)


def _walk(cpm, blk, varmap, boolmap, enforce=()):
    for c in blk.component_data_objects(Constraint, active=True, descend_into=False):
        _emit_constraint(cpm, c, varmap, enforce)

    if LogicalConstraint is not None:
        for lc in blk.component_data_objects(
            LogicalConstraint, active=True, descend_into=False
        ):
            _emit_logical(cpm, lc, boolmap, varmap, enforce)

    for disj in blk.component_data_objects(Disjunction, active=True, descend_into=False):
        inds = []
        for d in disj.disjuncts:
            if not d.active:
                continue
            ind = cpm.NewBoolVar(d.name)
            inds.append(ind)
            _walk(cpm, d, varmap, boolmap, tuple(enforce) + (ind,))
        _emit_selection(cpm, inds, bool(getattr(disj, "xor", True)), enforce)

    for sub in blk.component_data_objects(Block, active=True, descend_into=False):
        _walk(cpm, sub, varmap, boolmap, enforce)


def _emit_selection(cpm, inds, xor, enforce):
    if not inds:
        return
    lits = list(enforce)
    if not lits:
        if xor:
            cpm.AddExactlyOne(inds)
        else:
            cpm.AddBoolOr(inds)
    else:
        aux = cpm.NewBoolVar("_pc_sel")
        cpm.AddMinEquality(aux, lits)
        if xor:
            cpm.Add(sum(inds) == aux)
        else:
            cpm.Add(sum(inds) >= aux)


def build_cpsat_model(model):
    """Translate a Pyomo model into (CpModel, varmap, boolmap).

    varmap maps id(VarData) -> (pyomo_var, cp_var); boolmap maps
    id(BooleanVarData) -> cp_bool, both for solution load-back.
    """
    from ortools.sat.python import cp_model

    cpm = cp_model.CpModel()
    varmap = {}
    boolmap = {}

    for v in model.component_data_objects(Var, active=True, descend_into=True):
        if v.fixed or id(v) in varmap:
            continue
        if not (v.is_integer() or v.is_binary()):
            raise ValueError(
                f"pyomo-cp: variable '{v.name}' is continuous. The CP-SAT "
                f"backend is finite-domain; apply "
                f"TransformationFactory('cp.discretize') first."
            )
        lb, ub = v.bounds
        if lb is None or ub is None:
            raise ValueError(
                f"pyomo-cp: variable '{v.name}' must have finite bounds."
            )
        if v.is_binary():
            varmap[id(v)] = (v, cpm.NewBoolVar(v.name))
        else:
            varmap[id(v)] = (v, cpm.NewIntVar(
                _as_int(lb, "lower bound"), _as_int(ub, "upper bound"), v.name
            ))

    for bv in model.component_data_objects(BooleanVar, active=True, descend_into=True):
        _get_bool(cpm, bv, boolmap, varmap)

    _walk(cpm, model, varmap, boolmap)

    objs = [o for o in model.component_data_objects(Objective, active=True)]
    if len(objs) > 1:
        raise ValueError("pyomo-cp: multiple active objectives are not supported.")
    if objs:
        obj = objs[0]
        repn = generate_standard_repn(obj.expr)
        if not repn.is_linear():
            raise ValueError("pyomo-cp: the objective must be linear.")
        expr = _linear_expr(repn, varmap)
        if obj.sense == minimize:
            cpm.Minimize(expr)
        else:
            cpm.Maximize(expr)

    return cpm, varmap, boolmap


@SolverFactory.register(
    "cpsat", doc="CP-SAT (OR-Tools) backend for Pyomo (pyomo-cp)."
)
class CPSATSolver(OptSolver):
    """Solve a Pyomo model (integer, discretized, GDP, logical) with CP-SAT."""

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
        from ortools.sat.python import cp_model

        load_solutions = kwds.pop("load_solutions", True)
        time_limit = kwds.pop("time_limit", None)
        workers = kwds.pop("workers", None)
        seed = kwds.pop("seed", None)

        cpm, varmap, boolmap = build_cpsat_model(model)

        solver = cp_model.CpSolver()
        if time_limit is not None:
            solver.parameters.max_time_in_seconds = float(time_limit)
        if workers is not None:
            solver.parameters.num_search_workers = int(workers)
        if seed is not None:
            solver.parameters.random_seed = int(seed)
        status = solver.Solve(cpm)

        results = self._build_results(solver, status)
        if load_solutions and status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            for v, cpv in varmap.values():
                v.set_value(solver.Value(cpv))
            for bv in _iter_boolean_vars(model):
                if id(bv) in boolmap:
                    bv.set_value(bool(solver.Value(boolmap[id(bv)])))
        return results

    @staticmethod
    def _build_results(solver, status):
        from ortools.sat.python import cp_model

        tc = {
            cp_model.OPTIMAL: TerminationCondition.optimal,
            cp_model.FEASIBLE: TerminationCondition.feasible,
            cp_model.INFEASIBLE: TerminationCondition.infeasible,
            cp_model.MODEL_INVALID: TerminationCondition.error,
            cp_model.UNKNOWN: TerminationCondition.unknown,
        }.get(status, TerminationCondition.unknown)

        results = SolverResults()
        results.solver.name = "cpsat"
        results.solver.wallclock_time = solver.WallTime()
        results.solver.termination_condition = tc
        results.solver.status = (
            SolverStatus.ok
            if tc in (TerminationCondition.optimal, TerminationCondition.feasible)
            else SolverStatus.warning
        )
        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            try:
                results.problem.upper_bound = solver.ObjectiveValue()
                results.problem.lower_bound = solver.BestObjectiveBound()
            except Exception:  # noqa: BLE001
                pass
        return results


def _iter_boolean_vars(model):
    return model.component_data_objects(BooleanVar, active=True, descend_into=True)
