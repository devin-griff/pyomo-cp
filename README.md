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

## Planned usage

```python
import pyomo.environ as pyo
import pyomo_cp  # registers TransformationFactory('cp.discretize') and SolverFactory('cpsat')

m = pyo.ConcreteModel()
# ... build a pyomo.gdp model ...

pyo.TransformationFactory("cp.discretize").apply_to(m, step=1.0)  # explicit
results = pyo.SolverFactory("cpsat").solve(m, time_limit=60)
```

## License

Apache License 2.0. See [LICENSE](LICENSE).
