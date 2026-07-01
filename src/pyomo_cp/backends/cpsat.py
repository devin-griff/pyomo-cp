"""CP-SAT backend.

Translates a Pyomo model (variables, linear constraints, ``pyomo.gdp``
disjunctions, and a linear objective) into an OR-Tools CP-SAT model, solves it,
and loads the solution back onto the Pyomo variables.

Disjunctions map to native reified constraints: each disjunct gets an indicator
BoolVar, its constraints are enforced only when the indicator is true
(``OnlyEnforceIf``), and the disjunction adds an exactly-one / at-least-one
selection. Nesting is handled by accumulating the enclosing indicators.

Continuous variables are rejected with a pointer to
``TransformationFactory('cp.discretize')`` rather than discretized silently.

Registered as ``SolverFactory('cpsat')``. ``ortools`` is imported lazily so the
package imports without it; install ``pyomo-cp[cpsat]`` to use this backend.
"""
from pyomo.core import Block, Constraint, Objective, Var, minimize, value
from pyomo.gdp import Disjunction

try:  # LogicalConstraint location has moved across Pyomo versions
    from pyomo.core import LogicalConstraint
except ImportError:  # pragma: no cover
    LogicalConstraint = None
from pyomo.opt import SolverFactory, SolverResults, SolverStatus, TerminationCondition
from pyomo.opt.base.solvers import OptSolver
from pyomo.repn import generate_standard_repn

_INT_TOL = 1e-9


def _as_int(x, what="value"):
    """Return x as a Python int, or raise if it isn't (near-)integral."""
    xi = round(x)
    if abs(x - xi) > _INT_TOL:
        raise ValueError(
            f"pyomo-cp: the CP-SAT backend needs an integer {what}, got {x}. "
            f"Scale the data or apply TransformationFactory('cp.discretize')."
        )
    return int(xi)


def _linear_expr(repn, varmap):
    """CP-SAT linear expression from a Pyomo standard repn. varmap is keyed by
    id(var) -> (pyomo_var, cp_var) since Pyomo vars aren't hashable."""
    expr = _as_int(repn.constant, "constant")
    for coef, var in zip(repn.linear_coefs, repn.linear_vars):
        if id(var) in varmap:
            expr = expr + _as_int(coef, "coefficient") * varmap[id(var)][1]
        else:  # fixed var not folded into the constant: treat as a constant
            expr = expr + _as_int(coef * value(var), "coefficient")
    return expr


def _emit_constraint(cpm, c, varmap, enforce):
    """Add a linear constraint, reified on the enforcing indicators if any."""
    repn = generate_standard_repn(c.body)
    if not repn.is_linear():
        raise ValueError(
            f"pyomo-cp: constraint '{c.name}' is nonlinear; the CP-SAT backend "
            f"supports linear constraints only."
        )
    expr = _linear_expr(repn, varmap)
    lits = list(enforce)

    def add(relation):
        ct = cpm.Add(relation)
        if lits:
            ct.OnlyEnforceIf(lits)

    if c.equality:
        add(expr == _as_int(value(c.lower), "rhs"))
    else:
        if c.lower is not None:
            add(expr >= _as_int(value(c.lower), "lower bound"))
        if c.upper is not None:
            add(expr <= _as_int(value(c.upper), "upper bound"))


def _emit_selection(cpm, inds, xor, enforce):
    """Add the disjunction's selection: exactly-one (xor) or at-least-one, made
    conditional on the enclosing indicators when nested."""
    if not inds:
        return
    lits = list(enforce)
    if not lits:
        if xor:
            cpm.AddExactlyOne(inds)
        else:
            cpm.AddBoolOr(inds)
    else:
        # aux = AND(enforce) via min of booleans; then tie the count to it.
        aux = cpm.NewBoolVar("_pc_sel")
        cpm.AddMinEquality(aux, lits)
        if xor:
            cpm.Add(sum(inds) == aux)
        else:
            cpm.Add(sum(inds) >= aux)


def _walk(cpm, blk, varmap, enforce=()):
    """Recursively translate a block/disjunct's constraints and disjunctions,
    carrying the enclosing indicator literals in `enforce`."""
    for c in blk.component_data_objects(Constraint, active=True, descend_into=False):
        _emit_constraint(cpm, c, varmap, enforce)

    if LogicalConstraint is not None:
        for lc in blk.component_data_objects(
            LogicalConstraint, active=True, descend_into=False
        ):
            raise NotImplementedError(
                f"pyomo-cp: logical constraint '{lc.name}' is not supported yet "
                f"(Phase 2b). Express it via disjunctions for now."
            )

    for disj in blk.component_data_objects(Disjunction, active=True, descend_into=False):
        inds = []
        for d in disj.disjuncts:
            if not d.active:
                continue
            ind = cpm.NewBoolVar(d.name)
            inds.append(ind)
            _walk(cpm, d, varmap, tuple(enforce) + (ind,))
        _emit_selection(cpm, inds, bool(getattr(disj, "xor", True)), enforce)

    # regular sub-blocks (Disjuncts have ctype Disjunct, not Block, so excluded)
    for sub in blk.component_data_objects(Block, active=True, descend_into=False):
        _walk(cpm, sub, varmap, enforce)


def build_cpsat_model(model):
    """Translate a Pyomo model into (CpModel, varmap).

    varmap maps id(VarData) -> (pyomo_var, cp_var) for solution load-back.
    Raises ValueError on continuous/unbounded variables, nonlinear expressions,
    or non-integer data.
    """
    from ortools.sat.python import cp_model

    cpm = cp_model.CpModel()
    varmap = {}

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

    _walk(cpm, model, varmap)

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

    return cpm, varmap


@SolverFactory.register(
    "cpsat", doc="CP-SAT (OR-Tools) backend for Pyomo (pyomo-cp)."
)
class CPSATSolver(OptSolver):
    """Solve a Pyomo model (integer, discretized, and GDP) with OR-Tools CP-SAT."""

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

        cpm, varmap = build_cpsat_model(model)

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
            except Exception:
                pass
        return results
