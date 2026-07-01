"""Phase 3: cp.discretize + cpsat end-to-end on a continuous GDP model."""
import pytest
import pyomo.environ as pyo
from pyomo.gdp import Disjunct, Disjunction

import pyomo_cp  # noqa: F401

pytest.importorskip("ortools")


def _mini_layout():
    # Two 1D boxes of length 2 and 3 on [0, 10] that must not overlap; minimize
    # the extent L. Continuous positions. Best packing is adjacent -> L = 5.
    m = pyo.ConcreteModel()
    m.x1 = pyo.Var(bounds=(0, 10))  # continuous (Reals)
    m.x2 = pyo.Var(bounds=(0, 10))
    m.L = pyo.Var(bounds=(0, 10))
    m.e1 = pyo.Constraint(expr=m.L >= m.x1 + 2)
    m.e2 = pyo.Constraint(expr=m.L >= m.x2 + 3)
    m.d1 = Disjunct()
    m.d1.c = pyo.Constraint(expr=m.x1 + 2 <= m.x2)  # box1 left of box2
    m.d2 = Disjunct()
    m.d2.c = pyo.Constraint(expr=m.x2 + 3 <= m.x1)  # box2 left of box1
    m.no = Disjunction(expr=[m.d1, m.d2])
    m.obj = pyo.Objective(expr=m.L)  # minimize
    return m


def test_continuous_rejected_before_discretize():
    m = _mini_layout()
    with pytest.raises(ValueError, match="discretize"):
        pyo.SolverFactory("cpsat").solve(m)


def test_discretize_then_solve():
    m = _mini_layout()
    pyo.TransformationFactory("cp.discretize").apply_to(m)
    # after discretization the position variables are integer-domain
    assert m.x1.is_integer() and m.L.is_integer()
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m.obj) == 5
    # a valid non-overlapping packing at the optimum
    x1, x2 = pyo.value(m.x1), pyo.value(m.x2)
    assert (x1 + 2 <= x2) or (x2 + 3 <= x1)


def test_step_not_one_raises():
    m = _mini_layout()
    with pytest.raises(NotImplementedError):
        pyo.TransformationFactory("cp.discretize").apply_to(m, step=0.5)


def test_leaves_integer_vars_alone():
    m = pyo.ConcreteModel()
    m.n = pyo.Var(bounds=(0, 5), domain=pyo.Integers)
    m.obj = pyo.Objective(expr=m.n, sense=pyo.maximize)
    pyo.TransformationFactory("cp.discretize").apply_to(m)
    assert m.n.is_integer()
    res = pyo.SolverFactory("cpsat").solve(m)
    assert pyo.value(m.obj) == 5
