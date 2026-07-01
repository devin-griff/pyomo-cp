# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-07-01

First release (alpha).

### Added
- `cp.discretize` transformation (`TransformationFactory('cp.discretize')`):
  discretize bounded continuous variables onto an integer grid — the unit grid
  in place, or a non-unit `step` via `x = lb + step * x_int` substitution, with
  automatic descaling of the solution. Descends into `pyomo.gdp` disjuncts so
  constraints inside disjunctions are discretized too.
- Transform-time feasibility guard: a linear equality with no solution on the
  chosen grid (e.g. a dimension pinned to an odd value on an even-step grid, or
  `x + y == 5` on an even grid) is rejected at the transformation with a clear
  message, rather than surfacing as a puzzling `infeasible` at solve time.
- `cpsat` solver backend (`SolverFactory('cpsat')`): translate integer models,
  `pyomo.gdp` disjunctions (reified natively, no big-M), and logical constraints
  / Boolean variables to OR-Tools CP-SAT, solve, and load the solution back onto
  the Pyomo variables. Fractional constraint and objective coefficients are
  scaled to integers automatically.
- Solver options via friendly aliases (`time_limit`, `workers`, `seed`, `gap`)
  or raw CP-SAT parameter names through `options={...}`. `tee=True` streams the
  CP-SAT search log through Python's stdout (visible in notebooks and when stdout
  is redirected).

[Unreleased]: https://github.com/devin-griff/pyomo-cp/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/devin-griff/pyomo-cp/releases/tag/v0.1.0
