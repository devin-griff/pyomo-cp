# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
"""pyomo-cp: constraint-programming backends for Pyomo.

Importing this package registers the pyomo-cp plugins:

  * ``TransformationFactory('cp.discretize')`` — discretize bounded continuous
    variables onto integer grids (an explicit, solver-agnostic step).
  * ``SolverFactory('cpsat')`` — CP-SAT (OR-Tools) backend.

Scope for now is a *backend framework*: it translates what Pyomo can already
express (variables, linear constraints, logical constraints, and
``pyomo.gdp`` disjunctions) to CP solvers. Global-constraint modelling is a
possible future direction, not yet supported.
"""
from importlib.metadata import PackageNotFoundError, version

# Import for side effects: these modules register the plugins on import.
from . import discretize as _discretize  # noqa: F401
from .backends import cpsat as _cpsat  # noqa: F401

try:
    __version__ = version("pyomo-cp")
except PackageNotFoundError:  # not installed (e.g. running from a source tree)
    __version__ = "0.0.0"

__all__ = ["__version__"]
