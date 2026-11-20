"""Pyomo models written as .nl (with .col/.row), .lp and .mps by Pyomo's own
writers, to validate rusolve's .nl reader against its LP/MPS readers.

usage: python make_nl.py OUTDIR [features|big|all]
  features: ranges, equalities, one-sided rows, a constant in a constraint
            body, free / upper-only / fixed / negative-lower variables,
            general integers and binaries, a min objective with a constant.
  big:      10,000 variables x 2,000 rows, seeded random sparse MILP
            (ranged / >= / = rows, 200 integers), feasible by construction.
"""
import os, random, sys
import pyomo.environ as pe

out = sys.argv[1]
which = sys.argv[2] if len(sys.argv) > 2 else "all"
os.makedirs(out, exist_ok=True)

def write(m, stem):
    for ext in ("nl", "lp", "mps"):
        m.write(os.path.join(out, f"{stem}.{ext}"), io_options={"symbolic_solver_labels": True})

def features():
    m = pe.ConcreteModel()
    m.i = pe.Var(domain=pe.Integers, bounds=(-5, 5))
    m.b = pe.Var(domain=pe.Binary)
    m.f = pe.Var()                       # free
    m.u = pe.Var(bounds=(None, 4))       # upper only
    m.w = pe.Var(bounds=(1.5, 1.5))      # fixed
    m.p = pe.Var(domain=pe.NonNegativeReals)
    m.q = pe.Var(bounds=(-3, 2))
    m.obj = pe.Objective(expr=2*m.i - 3*m.b + m.f - m.u + 0.25*m.w + m.p + 4*m.q - 3.5)
    m.eq = pe.Constraint(expr=m.i + m.b + m.q == 2)
    m.rng = pe.Constraint(expr=pe.inequality(-1, m.f - m.p + m.u, 6))
    m.le = pe.Constraint(expr=m.f + 2*m.q <= 3)
    m.ge = pe.Constraint(expr=m.u - m.i >= -7)
    m.cst = pe.Constraint(expr=m.p + m.w + 10 >= 12)   # constant in the body
    m.dup = pe.Constraint(expr=m.i + m.i - m.b <= 9)   # repeated variable
    m.box = pe.Constraint(expr=pe.inequality(-20, m.f + m.u + m.p, 20))
    write(m, "features")

def big():
    rnd = random.Random(20261118)
    n, rows, nint = 10000, 2000, 200
    m = pe.ConcreteModel()
    m.J = pe.RangeSet(0, n - 1)
    m.x = pe.Var(m.J, bounds=(0, 10))
    for j in range(nint):
        m.x[j].domain = pe.Integers
        m.x[j].setlb(0); m.x[j].setub(10)
    x0 = [rnd.uniform(0, 10) if j >= nint else rnd.randint(0, 10) for j in range(n)]
    m.c = pe.ConstraintList()
    for i in range(rows):
        cols = rnd.sample(range(n), 8)
        coef = [rnd.choice([-3, -2, -1, 1, 2, 3, 4]) for _ in cols]
        act = sum(a * x0[j] for a, j in zip(coef, cols))
        e = sum(a * m.x[j] for a, j in zip(coef, cols))
        kind = i % 3
        if kind == 0:
            m.c.add(pe.inequality(act - 5, e, act + 5))
        elif kind == 1:
            m.c.add(e >= act - 3)
        else:
            m.c.add(e <= act + 3)
    m.obj = pe.Objective(expr=sum(rnd.uniform(-1, 1) * m.x[j] for j in range(n)) + 12.5)
    write(m, "big")

if which in ("features", "all"): features()
if which in ("big", "all"): big()
