* rusolve highs-oracle synthetic LP for the G1 ray fixtures (g1/rays):
* min -X0 - X1 + X2 s.t. R0: X0 - X1 <= 1, R1: X0 - X2 <= 2, R2: X1 + X2
* >= 1, X0, X1, X2 >= 0, X2 <= 4. Unbounded along X0 = X1 -> inf, so a
* primal ray exists.
NAME          RAYUNB
ROWS
 N  COST
 L  R0
 L  R1
 G  R2
COLUMNS
    X0        COST        -1.0         R0           1.0
    X0        R1           1.0
    X1        COST        -1.0         R0          -1.0
    X1        R2           1.0
    X2        COST         1.0         R1          -1.0
    X2        R2           1.0
RHS
    RHS       R0           1.0         R1           2.0
    RHS       R2           1.0
BOUNDS
 UP BND       X2           4.0
ENDATA
