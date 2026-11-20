* rusolve highs-oracle synthetic LP for the G1 ray fixtures (g1/rays):
* max X0 + 2 X1 - X2 + X3 s.t. E0: X0 + X1 + X2 = 4, R1: 6 <= X0 + X1 + X3
* <= 8 (RANGES), L4: X0 + X1 <= 1, X0, X1 >= 0, X2 in (-inf, 2], X3 in
* [0, 1]. E0 with L4 forces X2 >= 3, and R1 needs X0 + X1 + X3 >= 6:
* infeasible maximisation with an equality, a ranged row and a column with
* no lower bound, so the dual is unbounded.
NAME          RAYDUB
OBJSENSE
    MAX
ROWS
 N  COST
 E  E0
 G  R1
 L  L4
COLUMNS
    X0        COST         1.0         E0           1.0
    X0        R1           1.0         L4           1.0
    X1        COST         2.0         E0           1.0
    X1        R1           1.0         L4           1.0
    X2        COST        -1.0         E0           1.0
    X3        COST         1.0         R1           1.0
RHS
    RHS       E0           4.0         R1           6.0
    RHS       L4           1.0
RANGES
    RNG       R1           2.0
BOUNDS
 MI BND       X2
 UP BND       X2           2.0
 UP BND       X3           1.0
ENDATA
