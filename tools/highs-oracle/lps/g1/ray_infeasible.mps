* rusolve highs-oracle synthetic LP for the G1 ray fixtures (g1/rays):
* R0: X0 + X1 <= 1, R1: X0 + X1 + X2 >= 3, R2: X1 - X2 <= 0.5, X0, X1 in
* [0, 2], X2 in [0, 1]. Primal infeasible (R0 and R1 with X2 <= 1), so a
* dual ray exists.
NAME          RAYINF
ROWS
 N  COST
 L  R0
 G  R1
 L  R2
COLUMNS
    X0        COST         1.0         R0           1.0
    X0        R1           1.0
    X1        COST         2.0         R0           1.0
    X1        R1           1.0         R2           1.0
    X2        COST         1.0         R1           1.0
    X2        R2          -1.0
RHS
    RHS       R0           1.0         R1           3.0
    RHS       R2           0.5
BOUNDS
 UP BND       X0           2.0
 UP BND       X1           2.0
 UP BND       X2           1.0
ENDATA
