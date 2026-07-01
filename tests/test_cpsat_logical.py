"""Phase 2b: CP-SAT backend on logical constraints / Boolean variables."""
import pytest
import pyomo.environ as pyo
from pyomo.environ import BooleanVar, LogicalConstraint, atleast, lnot

import pyomo_cp  # noqa: F401

pytest.importorskip("ortools")


def _solve(m):
    res = pyo.SolverFactory("cpsat").solve(m)
    assert res.solver.termination_condition in (
        pyo.TerminationCondition.optimal,
        pyo.TerminationCondition.feasible,
    )
    return res


def test_implication():
    m = pyo.ConcreteModel()
    m.a = BooleanVar()
    m.b = BooleanVar()
    m.force_a = LogicalConstraint(expr=m.a)  # a true
    m.imp = LogicalConstraint(expr=m.a.implies(m.b))  # a -> b
    _solve(m)
    assert pyo.value(m.a) is True
    assert pyo.value(m.b) is True


def test_equivalence():
    m = pyo.ConcreteModel()
    m.p = BooleanVar()
    m.q = BooleanVar()
    m.fp = LogicalConstraint(expr=m.p)
    m.eq = LogicalConstraint(expr=m.p.equivalent_to(m.q))
    _solve(m)
    assert pyo.value(m.q) is True


def test_not_and_or():
    m = pyo.ConcreteModel()
    m.a = BooleanVar()
    m.b = BooleanVar()
    m.na = LogicalConstraint(expr=lnot(m.a))  # a false
    m.orab = LogicalConstraint(expr=m.a.lor(m.b))  # a or b -> b true
    _solve(m)
    assert pyo.value(m.a) is False
    assert pyo.value(m.b) is True


def test_atleast_top_level():
    m = pyo.ConcreteModel()
    m.x = BooleanVar([1, 2, 3])
    m.al = LogicalConstraint(expr=atleast(2, m.x[1], m.x[2], m.x[3]))
    _solve(m)
    assert sum(1 for i in (1, 2, 3) if pyo.value(m.x[i])) >= 2
