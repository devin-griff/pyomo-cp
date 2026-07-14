# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
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
    m.obj = pyo.Objective(expr=sum(val[i] * m.take[i] for i in m.I), sense=pyo.maximize)

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


def test_solver_options_passthrough():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 3), domain=pyo.Integers)
    m.obj = pyo.Objective(expr=m.x, sense=pyo.maximize)
    opt = pyo.SolverFactory("cpsat")
    # friendly aliases + a raw CP-SAT parameter
    res = opt.solve(m, time_limit=5, workers=2, options={"random_seed": 1})
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m.obj) == 3
    with pytest.raises(ValueError, match="unknown CP-SAT parameter"):
        opt.solve(m, options={"definitely_not_a_parameter": 1})


def test_tee_streams_through_python_stdout():
    # tee=True must route the CP-SAT log through sys.stdout (so it shows in
    # notebooks / respects redirection), not straight to the C-level fd 1.
    import contextlib
    import io

    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 5), domain=pyo.Integers)
    m.obj = pyo.Objective(expr=m.x, sense=pyo.maximize)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        pyo.SolverFactory("cpsat").solve(m, tee=True, time_limit=5)
    assert "CpSolverResponse summary" in buf.getvalue()

    quiet = io.StringIO()
    with contextlib.redirect_stdout(quiet):
        pyo.SolverFactory("cpsat").solve(m, time_limit=5)
    assert quiet.getvalue() == ""


def test_fractional_coefficients_scaled():
    # x + y/2 <= 2.5 over integer x,y in [0,5]; the backend scales by 2 to
    # 2x + y <= 5. Maximize x+y -> 5 (x=0, y=5).
    m = pyo.ConcreteModel()
    m.x = pyo.Var(bounds=(0, 5), domain=pyo.Integers)
    m.y = pyo.Var(bounds=(0, 5), domain=pyo.Integers)
    m.c = pyo.Constraint(expr=m.x + m.y / 2 <= 2.5)
    m.obj = pyo.Objective(expr=m.x + m.y, sense=pyo.maximize)
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m.obj) == 5


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
