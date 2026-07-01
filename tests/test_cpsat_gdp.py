"""Phase 2: CP-SAT backend on pyomo.gdp disjunctions, checked against known optima."""
import pytest
import pyomo.environ as pyo
from pyomo.gdp import Disjunct, Disjunction  # noqa: F401

import pyomo_cp  # noqa: F401

pytest.importorskip("ortools")


def test_xor_disjunction():
    # Either (x>=6, y<=1) or (y>=6, x<=1). Maximize x+y -> 11 via one branch.
    # If the disjunct constraints were enforced unconditionally the model would
    # be infeasible (x>=6 and x<=1), so a feasible optimum proves reification.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 10), domain=pyo.Integers)
    m.y = pyo.Var(bounds=(0, 10), domain=pyo.Integers)
    m.d1 = Disjunct()
    m.d1.a = pyo.Constraint(expr=m.x >= 6)
    m.d1.b = pyo.Constraint(expr=m.y <= 1)
    m.d2 = Disjunct()
    m.d2.a = pyo.Constraint(expr=m.y >= 6)
    m.d2.b = pyo.Constraint(expr=m.x <= 1)
    m.disj = Disjunction(expr=[m.d1, m.d2])
    m.obj = pyo.Objective(expr=m.x + m.y, sense=pyo.maximize)

    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m.obj) == 11
    x, y = pyo.value(m.x), pyo.value(m.y)
    assert (x >= 6 and y <= 1) or (y >= 6 and x <= 1)


def test_nested_disjunction():
    # outer1: x<=18 AND (x>=15 OR x<=2);  outer2: x==7.  Maximize x -> 18.
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 20), domain=pyo.Integers)
    m.o1 = Disjunct()
    m.o1.cap = pyo.Constraint(expr=m.x <= 18)
    m.o1.i1 = Disjunct()
    m.o1.i1.c = pyo.Constraint(expr=m.x >= 15)
    m.o1.i2 = Disjunct()
    m.o1.i2.c = pyo.Constraint(expr=m.x <= 2)
    m.o1.inner = Disjunction(expr=[m.o1.i1, m.o1.i2])
    m.o2 = Disjunct()
    m.o2.c = pyo.Constraint(expr=m.x == 7)
    m.outer = Disjunction(expr=[m.o1, m.o2])
    m.obj = pyo.Objective(expr=m.x, sense=pyo.maximize)

    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m.obj) == 18
