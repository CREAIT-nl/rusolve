* rusolve highs-oracle synthetic LP for the W1 IIS fixtures (w1/iis):
* R0: X0 + X1 >= 4 with 0 <= X0 <= 1, 0 <= X1 <= 2; R1: X2 - X0 <= 3 with
* X2 free. Infeasible; the IIS is R0 plus the upper bounds of X0 and X1
* (R1 and the free X2 are not in it).
NAME          IISBNDS
ROWS
 N  COST
 G  R0
 L  R1
COLUMNS
    X0        COST         1.0         R0           1.0
    X0        R1          -1.0
    X1        COST         1.0         R0           1.0
    X2        R1           1.0
RHS
    RHS       R0           4.0         R1           3.0
BOUNDS
 UP BND       X0           1.0
 UP BND       X1           2.0
 FR BND       X2
ENDATA
