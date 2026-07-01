"""Phase 1: CP-SAT backend on flat integer models, checked against known optima."""
import pytest
import pyomo.environ as pyo

import pyomo_cp  # noqa: F401  (registers the plugins)

ortools = pytest.importorskip("ortools")  # skip if OR-Tools isn't installed


def test_small_integer_lp():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 2), domain=pyo.Integers)
    m.y = pyo.Var(bounds=(0, 2), domain=pyo.Integers)
    m.c = pyo.Constraint(expr=m.x + m.y <= 3)
    m.obj = pyo.Objective(expr=m.x + m.y, sense=pyo.maximize)

    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m.obj) == 3
    assert pyo.value(m.x) + pyo.value(m.y) == 3


def test_binary_knapsack():
    # items (value, weight); capacity 5. Optimal value = 40 (items 2 and 3).
    val = {1: 10, 2: 25, 3: 15}
    wt = {1: 2, 2: 3, 3: 2}
    m = pyo.ConcreteModel()
    m.I = pyo.Set(initialize=[1, 2, 3])
    m.take = pyo.Var(m.I, domain=pyo.Binary)
    m.cap = pyo.Constraint(expr=sum(wt[i] * m.take[i] for i in m.I) <= 5)
    m.obj = pyo.Objective(
        expr=sum(val[i] * m.take[i] for i in m.I), sense=pyo.maximize
    )

    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m.obj) == 40
    assert pyo.value(m.take[2]) == 1 and pyo.value(m.take[3]) == 1


def test_equality_and_minimize():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 10), domain=pyo.Integers)
    m.y = pyo.Var(bounds=(0, 10), domain=pyo.Integers)
    m.eq = pyo.Constraint(expr=m.x + 2 * m.y == 8)
    m.obj = pyo.Objective(expr=m.x + m.y)  # minimize
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    # minimize x+y with x+2y=8, x,y>=0 -> y=4, x=0, obj=4
    assert pyo.value(m.obj) == 4


def test_continuous_variable_rejected():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 1))  # continuous (Reals)
    m.obj = pyo.Objective(expr=m.x, sense=pyo.maximize)
    with pytest.raises(ValueError, match="discretize"):
        pyo.SolverFactory("cpsat").solve(m)


def test_infeasible():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 2), domain=pyo.Integers)
    m.c1 = pyo.Constraint(expr=m.x >= 5)  # impossible given bounds
    m.obj = pyo.Objective(expr=m.x)
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.infeasible
