# pyomo-cp

Constraint-programming backends for [Pyomo](https://www.pyomo.org/): a
discretization transform plus CP solver interfaces, starting with
[CP-SAT](https://developers.google.com/optimization/cp/cp_solver) (OR-Tools).

Pyomo can build and solve disjunctive/logical models via `pyomo.gdp`, but it can
only *solve* them by reformulating to MILP/MINLP (big-M, hull) and calling a MIP
solver. `pyomo-cp` adds the missing path: take the same model and solve it with a
constraint-programming solver, where disjunctions map to native reified
constraints instead of being reformulated.

> **Status: alpha.** Integer models, `pyomo.gdp` disjunctions, logical
> constraints / Boolean variables, and discretization of continuous variables
> (any `step`) work end-to-end via `SolverFactory('cpsat')`. Global constraints
> are not yet supported. See [ROADMAP.md](ROADMAP.md).

## Scope

`pyomo-cp` is a **backend framework**, not a full CP modelling frontend. It
translates what Pyomo can already express:

- integer and (once discretized) continuous variables,
- linear constraints and a linear objective,
- logical constraints (`BooleanVar` / `LogicalConstraint`),
- `pyomo.gdp` disjunctions.

**Global constraints** (`alldifferent`, `no_overlap`, `element`, `cumulative`,
...) are the heart of CP's modelling power, and Pyomo has no vocabulary for
them. Adding one is roadmap, not part of the initial scope.

CP solvers are finite-domain, so continuous variables must be discretized first.
That is an **explicit** step (`TransformationFactory('cp.discretize')`), never
automatic, because the integrality assumption is a modelling decision that
changes the problem.

## Install

```bash
pip install "pyomo-cp[cpsat]"   # includes OR-Tools for the CP-SAT backend
```

## Usage

```python
import pyomo.environ as pyo
from pyomo.gdp import Disjunct, Disjunction
import pyomo_cp  # registers cp.discretize and the cpsat solver

# Two boxes (lengths 2 and 3) that must not overlap on a line; minimize extent.
m = pyo.ConcreteModel()
m.x1 = pyo.Var(bounds=(0, 10))          # continuous positions
m.x2 = pyo.Var(bounds=(0, 10))
m.L = pyo.Var(bounds=(0, 10))
m.e1 = pyo.Constraint(expr=m.L >= m.x1 + 2)
m.e2 = pyo.Constraint(expr=m.L >= m.x2 + 3)
m.d1 = Disjunct(); m.d1.c = pyo.Constraint(expr=m.x1 + 2 <= m.x2)  # box1 left of box2
m.d2 = Disjunct(); m.d2.c = pyo.Constraint(expr=m.x2 + 3 <= m.x1)  # box2 left of box1
m.no = Disjunction(expr=[m.d1, m.d2])
m.obj = pyo.Objective(expr=m.L)

pyo.TransformationFactory("cp.discretize").apply_to(m)      # explicit; unit grid
res = pyo.SolverFactory("cpsat").solve(m, time_limit=10)
print(res.solver.termination_condition, pyo.value(m.obj))  # optimal 5.0
```

Solver options are passed as friendly aliases (`time_limit`, `workers`, `seed`,
`gap`) or as raw CP-SAT parameter names via `options={...}`, e.g.
`solve(m, workers=8, options={"log_search_progress": True})`.

## License

Apache License 2.0. See [LICENSE](LICENSE).
