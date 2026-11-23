* rusolve highs-oracle synthetic LP for the W1 IIS fixtures (w1/iis):
* R0: X0 + X1 <= 1, R1: X0 - X1 >= 2, X >= 0. Infeasible; the IIS is
* rows {R0, R1} plus the lower bound of X1 (X0 >= 0 is not needed).
* Kept out of lps/*.mps so the IPX generation loop does not pick it up.
NAME          IISROWS
ROWS
 N  COST
 L  R0
 G  R1
COLUMNS
    X0        R0           1.0         R1           1.0
    X1        R0           1.0         R1          -1.0
RHS
    RHS       R0           1.0         R1           2.0
ENDATA
