# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.2.0] - 2026-07-14

### Changed

- Relicensed from Apache-2.0 to BSD-3-Clause, matching the Pyomo and
  scientific-Python ecosystem this builds on.
- Minimum Python is now 3.10 (3.9 is end-of-life), and the minimum Pyomo is
  corrected to 6.8.0: the backend has always needed the public Data-class names
  (from 6.7.2) and NumPy 2.0 support (from 6.8.0), so the old `>=6.7` pin never
  actually worked.
- `cp.discretize` validates its `step` option through a `ConfigDict`, so an
  unknown keyword option now raises `ValueError` instead of being ignored.

### Internal

- Aligned with Pyomo's contribution conventions: NumPy-style docstrings on the
  public and private API, Black formatting with Pyomo's settings, per-file BSD
  license headers, and `attempt_import` for the optional OR-Tools dependency.
- CI gained Black and spell-check gates, coverage, a no-OR-Tools import check,
  and a minimum-dependency job, and it cuts a GitHub Release from the changelog
  on each version tag.

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

[Unreleased]: https://github.com/devin-griff/pyomo-cp/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/devin-griff/pyomo-cp/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/devin-griff/pyomo-cp/releases/tag/v0.1.0
