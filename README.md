# pyomo-cp

[![PyPI](https://img.shields.io/pypi/v/pyomo-cp.svg)](https://pypi.org/project/pyomo-cp/)
[![Python versions](https://img.shields.io/pypi/pyversions/pyomo-cp.svg)](https://pypi.org/project/pyomo-cp/)
[![CI](https://github.com/devin-griff/pyomo-cp/actions/workflows/ci.yml/badge.svg)](https://github.com/devin-griff/pyomo-cp/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

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
> are not yet supported.

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

## Relationship to `pyomo.contrib.cp`

Pyomo ships an in-tree constraint-programming module, `pyomo.contrib.cp`. It is a
different tool for a different job, and the two are complementary:

| | `pyomo.contrib.cp` | this package |
|---|---|---|
| Role | CP **frontend** | CP **backend** |
| Built for | scheduling: `IntervalVar`, sequence vars, `Pulse` / step functions, precedence | solving models Pyomo already expresses — no new constructs |
| Input | a CP model written with the interval API | an existing `pyomo.gdp` / integer / logical model |
| Continuous variables | none (finite-domain only) | explicit `cp.discretize` onto a grid |
| Solver | IBM CP Optimizer (commercial, via `docplex`) | CP-SAT / OR-Tools (open-source) |

In one line: `pyomo.contrib.cp` is *write an interval/scheduling model and solve
it with CP Optimizer*; this package is *take a disjunctive/integer model you'd
otherwise reformulate to MILP and solve it (optionally discretized) with CP-SAT*.
They overlap only on logical constraints; the paradigm (interval scheduling vs
disjunctive/geometric) and the solver (commercial vs open) otherwise differ.

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

See [examples/facility_layout.ipynb](examples/facility_layout.ipynb) for a
complete worked example: a 17-block plant-layout GDP model solved to proven
optimality by CP-SAT via `cp.discretize`.

Solver options are passed as friendly aliases (`time_limit`, `workers`, `seed`,
`gap`) or as raw CP-SAT parameter names via `options={...}`, e.g.
`solve(m, workers=8, options={"log_search_progress": True})`.

## License

Apache License 2.0. See [LICENSE](LICENSE).
