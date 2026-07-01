# pyomo-cp roadmap

The v0.1 goal: a user writes a `pyomo.gdp` model, applies `cp.discretize`, calls
`SolverFactory('cpsat').solve(m)`, and gets the correct optimum, validated
against an independent hand-written CP-SAT model and against Gurobi/MILP.

## Guiding principles

- **Backend framework, not a CP frontend.** Translate what Pyomo can already
  express; global-constraint modelling is out of scope for v0.
- **Oracle-tested.** Every translation is validated by solving the same model two
  independent ways and asserting identical optima.
- **Don't over-abstract.** Build CP-SAT as a direct emitter with a clean seam.
  Introduce a solver-agnostic IR only when a second backend forces it.
- **Discretization is explicit and general.** A first-class, solver-agnostic
  transform the user applies knowingly; backends never discretize silently.

## Scope contract (v0.1)

Supported: integer/continuous-bounded `Var`; linear `Constraint`; linear
`Objective`; `pyomo.gdp` `Disjunction`/`Disjunct`; `LogicalConstraint`/
`BooleanVar`. Requires all variables bounded and data integer-scalable (exact) or
discretized (approximate, declared `step`). Must error clearly on: nonlinear
expressions, unbounded variables, global constraints.

## Phase 0 — Scaffolding (done)

Repo, `src/pyomo_cp/` layout, `pyproject.toml` (Apache-2.0), stub
`cp.discretize` transform and `cpsat` solver, CI skeleton.
**Acceptance:** `pip install -e .` works; importing `pyomo_cp` registers
`SolverFactory('cpsat')` and `TransformationFactory('cp.discretize')`.

## Phase 1 — CP-SAT backend for already-integer models

Model walker over `Var`/`Constraint`/`Objective`; linear extraction via
`generate_standard_repn` (reject nonlinear with a clear error); emit
`NewIntVar`/`NewBoolVar`, `m.Add(...)`, `Minimize`/`Maximize`; solve wrapper
(time limit, workers, seed), status mapping, solution load-back, Pyomo results
object.
**Acceptance:** a small pure-integer model solves via `cpsat` and matches Gurobi.

## Phase 2 — Disjunctions and logic

`Disjunction` -> one indicator `BoolVar` per `Disjunct`, disjunct constraints
under `OnlyEnforceIf`, `AddExactlyOne`/`AddAtMostOne` per semantics (nested
handled recursively). `LogicalConstraint`/`BooleanVar` -> CP-SAT boolean
constraints.
**Acceptance:** a small disjunctive integer model matches its Gurobi (bigm/hull)
optimum.

## Phase 3 — Discretization transform

`cp.discretize`: bounded continuous `Var` -> integer `x_int` with
`x = lb + step*x_int`; substitute; scale coefficients to integers (exact for
integer data, warn/round otherwise); **preserve disjunctions**; descale on
load-back; per-variable `step`.
**Acceptance (headline milestone):** a `pyomo.gdp` layout model, discretized and
solved with `cpsat`, reproduces the hand-written CP-SAT optima on the benchmark
instances and equals the Gurobi optima. This is the "one GDP model, two solvers"
demonstration.

## Phase 4 — Hardening + release v0.1

Clear errors for unsupported constructs; solver-option passthrough; README
example; docs on scope and the discretization contract; green CI matrix; PyPI
trusted-publishing on tag; `v0.1.0`.

## Cross-cutting

- **Testing:** oracle tests (vs hand CP-SAT and vs Gurobi) plus property tests
  (random small models solved both ways with equal optima). Test solution
  *values*, not just objectives.
- **Version drift:** CI-pin ortools/pyomo; the interface is what breaks on their
  releases.

## Post-v0 (future)

1. **MiniZinc/FlatZinc emitter** — reaches Chuffed/Gecode/etc.; a very different
   target than CP-SAT's API, so it's what forces and validates a real IR. Factor
   the IR here, not before.
2. **Global-constraint vocabulary** — modelling components (`AllDifferent`,
   `NoOverlap`, `Element`, `Cumulative`) mapped to CP-SAT and MiniZinc globals;
   includes recognizing a GDP no-overlap pattern to emit `AddNoOverlap2D`. This
   is the "backend framework -> CP frontend" leap and effectively its own
   project.
3. **CP Optimizer backend** (docplex), once the IR exists.
4. **Ergonomics:** warm starts, search hints, portfolio/worker config.

## Open decisions

- Introduce the solver-agnostic IR at Phase 5 (MiniZinc), not before.
- License: Apache-2.0 (chosen: permissive, patent grant, matches OR-Tools).
