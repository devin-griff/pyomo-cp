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


def test_general_grid_step_half():
    # minimize x, continuous in [0,5], x >= 2.3, on a step-0.5 grid.
    # x = 0.5*xi, xi in [0,10]; x>=2.3 -> xi>=5 -> x = 2.5.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 5))
    m.c = pyo.Constraint(expr=m.x >= 2.3)
    m.obj = pyo.Objective(expr=m.x)  # minimize
    pyo.TransformationFactory("cp.discretize").apply_to(m, step=0.5)
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert abs(pyo.value(m.x) - 2.5) < 1e-9
    assert abs(pyo.value(m.obj) - 2.5) < 1e-9
    assert abs(res.problem.upper_bound - 2.5) < 1e-9


def test_offgrid_pin_raises_at_transform():
    # A variable pinned (single-var equality) to a value off the grid must fail
    # at the transformation, not surface as a puzzling infeasible at solve time.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 10))
    m.fix = pyo.Constraint(expr=m.x == 3)  # 3 is not on the even (step-2) grid
    m.obj = pyo.Objective(expr=m.x)
    with pytest.raises(ValueError, match="not on its discretization grid"):
        pyo.TransformationFactory("cp.discretize").apply_to(m, step=2)


def test_offgrid_pin_unit_grid_noninteger():
    # Same guard on the unit grid: a non-integer pin is caught at transform.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 10))
    m.fix = pyo.Constraint(expr=m.x == 2.5)
    m.obj = pyo.Objective(expr=m.x)
    with pytest.raises(ValueError, match="not on its discretization grid"):
        pyo.TransformationFactory("cp.discretize").apply_to(m)  # step=1


def test_ongrid_pin_ok():
    # A pin that lands on the grid transforms and solves fine.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 10))
    m.fix = pyo.Constraint(expr=m.x == 4)  # on the even grid
    m.obj = pyo.Objective(expr=m.x)
    pyo.TransformationFactory("cp.discretize").apply_to(m, step=2)
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert abs(pyo.value(m.x) - 4) < 1e-9


def test_multivar_sum_equality_unreachable_raises():
    # Each variable fits the even grid, but the sum can't be odd: x+y==5 has no
    # grid solution. Caught at transform by the gcd/divisibility test.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 10))
    m.y = pyo.Var(bounds=(0, 10))
    m.c = pyo.Constraint(expr=m.x + m.y == 5)
    m.obj = pyo.Objective(expr=m.x)
    with pytest.raises(ValueError, match="no solution on the discretization grid"):
        pyo.TransformationFactory("cp.discretize").apply_to(m, step=2)


def test_multivar_sum_equality_reachable_ok():
    # x+y==4 is reachable on the even grid; transforms and solves.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 10))
    m.y = pyo.Var(bounds=(0, 10))
    m.c = pyo.Constraint(expr=m.x + m.y == 4)
    m.obj = pyo.Objective(expr=m.x, sense=pyo.maximize)
    pyo.TransformationFactory("cp.discretize").apply_to(m, step=2)
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert abs(pyo.value(m.x) + pyo.value(m.y) - 4) < 1e-9


def test_invalid_step_raises():
    m = _mini_layout()
    with pytest.raises(ValueError):
        pyo.TransformationFactory("cp.discretize").apply_to(m, step=0)


def test_leaves_integer_vars_alone():
    m = pyo.ConcreteModel()
    m.n = pyo.Var(bounds=(0, 5), domain=pyo.Integers)
    m.obj = pyo.Objective(expr=m.n, sense=pyo.maximize)
    pyo.TransformationFactory("cp.discretize").apply_to(m)
    assert m.n.is_integer()
    res = pyo.SolverFactory("cpsat").solve(m)
    assert pyo.value(m.obj) == 5
