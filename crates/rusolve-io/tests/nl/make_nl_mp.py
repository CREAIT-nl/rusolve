"""The tests/nl/features model written by AMPL's own NL writer (the AMPL MP
library's NLWriter2, github.com/ampl/mp, through its `nlwpy` Python
bindings, BSD-3, `pip install nlwpy`), as binary and as text .nl with
.col/.row name files. This is the real-producer check for binary .nl:
Pyomo writes text only.

usage: python make_nl_mp.py OUTDIR [MPS]
  without MPS: features_mp_binary.nl / features_mp_text.nl (+ .col/.row)
  with MPS:    OUTDIR/<stem>_mp_binary.nl from that MPS file read by highspy
               (for large-model checks; not a committed fixture)
"""
import os, sys
import nlwpy as n

INF = float("inf")
out = sys.argv[1]
os.makedirs(out, exist_ok=True)


def write(stem, cols, rows, obj, text):
    """cols: [(name, lb, ub, is_int)]; rows: [(name, lo, hi, {j: a})];
    obj: (sense, constant, [c_j]). The writer wants integers last."""
    m = n.NLModel(stem)
    m.SetCols([c[1] for c in cols], [c[2] for c in cols],
              [n.VarType.Integer if c[3] else n.VarType.Continuous for c in cols])
    m.SetColNames([c[0] for c in cols])
    starts, idx, val = [], [], []
    for r in rows:
        starts.append(len(idx))
        for j, a in sorted(r[3].items()):
            idx.append(j)
            val.append(float(a))
    m.SetRows([r[1] for r in rows], [r[2] for r in rows], n.MatrixFormat.Rowwise,
              starts, idx, val)
    m.SetRowNames([r[0] for r in rows])
    m.SetLinearObjective(obj[0], obj[1], obj[2])
    m.SetObjName("obj")
    s = n.NLSolver()
    o = n.MakeNLOptionsBasic_Default()
    o.n_text_mode_ = int(text)
    o.want_nl_comments_ = int(text)
    s.SetNLOptions(o)
    s.SetFileStub(os.path.join(out, stem))
    # Writes stem.nl/.col/.row, then runs the "solver"; /usr/bin/true writes
    # no .sol, so the reported error is expected and ignored.
    s.Solve(m, "/usr/bin/true", "")


def features():
    cols = [("f", -INF, INF, 0), ("u", -INF, 4, 0), ("w", 1.5, 1.5, 0),
            ("p", 0, INF, 0), ("q", -3, 2, 0), ("b", 0, 1, 1), ("i", -5, 5, 1)]
    rows = [("eq", 2, 2, {6: 1, 5: 1, 4: 1}), ("rng", -1, 6, {0: 1, 1: 1, 3: -1}),
            ("le", -INF, 3, {0: 1, 4: 2}), ("ge", -7, INF, {1: 1, 6: -1}),
            ("cst", 2, INF, {2: 1, 3: 1}), ("dup", -INF, 9, {6: 2, 5: -1}),
            ("box", -20, 20, {0: 1, 1: 1, 3: 1})]
    obj = (n.ObjSense.Minimize, -3.5, [1, -1, 0.25, 1, 4, -3, 2])
    write("features_mp_binary", cols, rows, obj, False)
    write("features_mp_text", cols, rows, obj, True)


def from_mps(path):
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.readModel(path)
    lp = h.getLp()
    nc, nr = lp.num_col_, lp.num_row_
    integ = list(lp.integrality_) or [highspy.HighsVarType.kContinuous] * nc
    is_int = [t == highspy.HighsVarType.kInteger for t in integ]
    order = [j for j in range(nc) if not is_int[j]] + [j for j in range(nc) if is_int[j]]
    pos = {j: k for k, j in enumerate(order)}
    cols = [(lp.col_names_[j], lp.col_lower_[j], lp.col_upper_[j], is_int[j]) for j in order]
    a = lp.a_matrix_
    rows = [(lp.row_names_[i], lp.row_lower_[i], lp.row_upper_[i], {}) for i in range(nr)]
    for j in range(nc):
        for k in range(a.start_[j], a.start_[j + 1]):
            rows[a.index_[k]][3][pos[j]] = a.value_[k]
    sense = n.ObjSense.Maximize if lp.sense_ == highspy.ObjSense.kMaximize else n.ObjSense.Minimize
    obj = (sense, lp.offset_, [lp.col_cost_[j] for j in order])
    stem = os.path.basename(path).split(".")[0] + "_mp_binary"
    write(stem, cols, rows, obj, False)


if len(sys.argv) > 2:
    from_mps(sys.argv[2])
else:
    features()
