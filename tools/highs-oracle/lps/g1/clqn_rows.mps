* G1 clique-normalise (HiGHS 1.16.0 HPresolveProbingUtil.cpp:16-176,
* normaliseCliqueRows): one row per branch, each on its own columns.
* Every column is binary except p0 (integer, [0,3]) and y (continuous, [0,1]).
NAME          clqn_rows
ROWS
 N  obj
 L  c00_unequal
 L  c01_setppc
 G  c02_geq_comp
 L  c03_fix_no_clique
 L  c04_complemented
 E  c05_equation
 E  c06_equation_unequal
 G  c07_ranged_setppc
 G  c08_ranged_equation
 G  c09_ranged_lhs
 L  c10_no_clique
 L  c11_fractional
 L  c12_dynamism
 L  c13_general_int
 L  c14_continuous
 L  c15_redundant
 L  c16_setppc_comp
 G  c17_cover
 G  c18_geq_clique
 L  c19_mixed_sign
 L  c20_fix_then_clique
 L  c21_singleton
 E  c22_equation_comp
 G  c23_ranged_both
 G  c24_ranged_rhs
 G  c25_geq_fix_clique
 E  c26_equation_scaled_comp
COLUMNS
    MARKER                 'MARKER'                 'INTORG'
    a0        obj       1
    a0        c00_unequal               2
    a1        c00_unequal               2
    a2        obj       1
    a2        c00_unequal               3
    b0        c01_setppc                1
    b1        obj       1
    b1        c01_setppc                1
    b2        c01_setppc                1
    d0        obj       1
    d0        c02_geq_comp              -1
    d1        c02_geq_comp              -1
    e0        obj       1
    e0        c03_fix_no_clique         3
    e1        c03_fix_no_clique         1
    e2        obj       1
    e2        c03_fix_no_clique         1
    f0        c04_complemented          3
    f1        obj       1
    f1        c04_complemented          -1
    f2        c04_complemented          -1
    g0        obj       1
    g0        c05_equation              2
    g1        c05_equation              2
    g2        obj       1
    g2        c05_equation              2
    h0        c06_equation_unequal      2
    h1        obj       1
    h1        c06_equation_unequal      2
    h2        c06_equation_unequal      3
    i0        obj       1
    i0        c07_ranged_setppc         1
    i1        c07_ranged_setppc         1
    j0        obj       1
    j0        c08_ranged_equation       2
    j1        c08_ranged_equation       2
    k0        obj       1
    k0        c09_ranged_lhs            -1
    k3        c09_ranged_lhs            1
    l0        obj       1
    l0        c10_no_clique             2
    l1        c10_no_clique             3
    l2        obj       1
    l2        c10_no_clique             4
    m0        c11_fractional            2.5
    m1        obj       1
    m1        c11_fractional            2.5
    n0        c12_dynamism              1
    n1        obj       1
    n1        c12_dynamism              200000
    o0        c13_general_int           1
    o1        obj       1
    o1        c13_general_int           1
    p0        c13_general_int           1
    q0        obj       1
    q0        c14_continuous            1
    q1        c14_continuous            1
    MARKER                 'MARKER'                 'INTEND'
    y         obj       1
    y         c14_continuous            1
    MARKER                 'MARKER'                 'INTORG'
    r0        c15_redundant             1
    r1        obj       1
    r1        c15_redundant             1
    s0        c16_setppc_comp           -1
    s1        obj       1
    s1        c16_setppc_comp           -1
    s2        c16_setppc_comp           -1
    t0        obj       1
    t0        c17_cover                 3
    t1        c17_cover                 3
    t2        obj       1
    t2        c17_cover                 3
    u0        c18_geq_clique            4
    u1        obj       1
    u1        c18_geq_clique            4
    u2        c18_geq_clique            4
    v0        obj       1
    v0        c19_mixed_sign            2
    v1        c19_mixed_sign            -3
    v2        obj       1
    v2        c19_mixed_sign            1
    w0        c20_fix_then_clique       5
    w1        obj       1
    w1        c20_fix_then_clique       3
    w2        c20_fix_then_clique       3
    x0        obj       1
    x0        c21_singleton             2
    z0        c22_equation_comp         -1
    z1        obj       1
    z1        c22_equation_comp         -1
    a3        c23_ranged_both           2
    a4        obj       1
    a4        c23_ranged_both           2
    a5        c23_ranged_both           3
    aa0       obj       1
    aa0       c24_ranged_rhs            1
    aa1       c24_ranged_rhs            -2
    aa2       obj       1
    aa2       c24_ranged_rhs            -2
    ab0       c25_geq_fix_clique        -4
    ab1       obj       1
    ab1       c25_geq_fix_clique        -2
    ab2       c25_geq_fix_clique        -2
    ac0       obj       1
    ac0       c26_equation_scaled_comp  -2
    ac1       c26_equation_scaled_comp  -2
    ac2       obj       1
    ac2       c26_equation_scaled_comp  2
    MARKER                 'MARKER'                 'INTEND'
RHS
    RHS       c00_unequal               3
    RHS       c01_setppc                1
    RHS       c02_geq_comp              -1
    RHS       c03_fix_no_clique         2
    RHS       c04_complemented          -1
    RHS       c05_equation              2
    RHS       c06_equation_unequal      3
    RHS       c07_ranged_setppc         0.5
    RHS       c08_ranged_equation       1
    RHS       c09_ranged_lhs            0
    RHS       c10_no_clique             5
    RHS       c11_fractional            3.7
    RHS       c12_dynamism              150000
    RHS       c13_general_int           1
    RHS       c14_continuous            1
    RHS       c15_redundant             5
    RHS       c16_setppc_comp           -2
    RHS       c17_cover                 3
    RHS       c18_geq_clique            8
    RHS       c19_mixed_sign            1
    RHS       c20_fix_then_clique       4
    RHS       c21_singleton             1
    RHS       c22_equation_comp         -1
    RHS       c23_ranged_both           -1
    RHS       c24_ranged_rhs            1
    RHS       c25_geq_fix_clique        -3
    RHS       c26_equation_scaled_comp  -2
RANGES
    RNG       c07_ranged_setppc         0.5
    RNG       c08_ranged_equation       1
    RNG       c09_ranged_lhs            1
    RNG       c23_ranged_both           4
BOUNDS
 BV BND       a0
 BV BND       a1
 BV BND       a2
 BV BND       b0
 BV BND       b1
 BV BND       b2
 BV BND       d0
 BV BND       d1
 BV BND       e0
 BV BND       e1
 BV BND       e2
 BV BND       f0
 BV BND       f1
 BV BND       f2
 BV BND       g0
 BV BND       g1
 BV BND       g2
 BV BND       h0
 BV BND       h1
 BV BND       h2
 BV BND       i0
 BV BND       i1
 BV BND       j0
 BV BND       j1
 BV BND       k0
 BV BND       k3
 BV BND       l0
 BV BND       l1
 BV BND       l2
 BV BND       m0
 BV BND       m1
 BV BND       n0
 BV BND       n1
 BV BND       o0
 BV BND       o1
 LO BND       p0        0
 UP BND       p0        3
 BV BND       q0
 BV BND       q1
 LO BND       y         0
 UP BND       y         1
 BV BND       r0
 BV BND       r1
 BV BND       s0
 BV BND       s1
 BV BND       s2
 BV BND       t0
 BV BND       t1
 BV BND       t2
 BV BND       u0
 BV BND       u1
 BV BND       u2
 BV BND       v0
 BV BND       v1
 BV BND       v2
 BV BND       w0
 BV BND       w1
 BV BND       w2
 BV BND       x0
 BV BND       z0
 BV BND       z1
 BV BND       a3
 BV BND       a4
 BV BND       a5
 BV BND       aa0
 BV BND       aa1
 BV BND       aa2
 BV BND       ab0
 BV BND       ab1
 BV BND       ab2
 BV BND       ac0
 BV BND       ac1
 BV BND       ac2
ENDATA
