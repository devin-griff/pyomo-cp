"""Phase 0 smoke tests: the package imports and registers its plugins."""
import pyomo_cp  # noqa: F401  (import registers the plugins)
from pyomo.core import TransformationFactory
from pyomo.opt import SolverFactory


def test_version_is_string():
    assert isinstance(pyomo_cp.__version__, str)


def test_discretize_transformation_registered():
    assert TransformationFactory("cp.discretize") is not None


def test_cpsat_solver_registered():
    assert SolverFactory("cpsat") is not None
