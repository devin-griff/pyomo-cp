"""CP solver backends for pyomo-cp.

CP-SAT (OR-Tools) is the first backend. A MiniZinc/FlatZinc emitter (reaching
Chuffed, Gecode, and others) and an IBM CP Optimizer backend are roadmap items;
see ROADMAP.md. The solver-agnostic intermediate representation those backends
would share is intentionally deferred until a second backend exists to validate
its shape.
"""
