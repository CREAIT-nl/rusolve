* rusolve highs-oracle synthetic LP for the IPX fixtures: free (x1), boxed
* (x2, x6 with negative lb), fixed (x3), upper-bounded only (x4, flipped by
* ipx::Model), shifted lower bound (x5), ranged row (R5), badly scaled row
* (R4, exercises equilibration).
NAME          IPXMIX
ROWS
 N  COST
 G  R1
 G  R2
 E  R3
 L  R4
 G  R5
COLUMNS
    X1        COST         1.0         R1           1.0
    X1        R3           1.0         R5           1.0
    X2        COST         2.0         R2           1.0
    X2        R3           1.0         R4        1000.0
    X3        COST        -1.0         R3           1.0
    X4        COST         0.5         R1          -1.0
    X4        R2           1.0
    X5        COST         3.0         R1           1.0
    X5        R4         0.001
    X6        COST        -1.0         R3           1.0
    X6        R5          -1.0
RHS
    RHS       R1          -1.0         R2          -3.0
    RHS       R3           6.0         R4        3000.0
    RHS       R5          -5.0
RANGES
    RNG       R5          13.0
BOUNDS
 FR BND       X1
 UP BND       X2           4.0
 FX BND       X3           2.0
 MI BND       X4
 UP BND       X4           3.0
 LO BND       X5           1.0
 LO BND       X6          -2.0
 UP BND       X6           5.0
ENDATA
