* rusolve highs-oracle synthetic LP for the G1 standard-form fixtures
* (g1/standard_form): max 2 X1 - 3 X2 + X3 + 4 X4 + X5 + 7 s.t.
* R1: X1 + X2 + X4 <= 10, R2: X1 - X3 + X5 >= -2, R3: X2 + X4 + X5 = 4,
* R4: -5 <= X1 - X4 + 2 X5 <= 8 (RANGES), R5: X1 + X2 + X3 + X4 <= 12;
* X1 free, X2 in [0, 5], X3 = 2, X4 in [-3, 6], X5 >= -1. The objective
* RHS -7 is the offset +7.
NAME          STDFORM
OBJSENSE
    MAX
ROWS
 N  COST
 L  R1
 G  R2
 E  R3
 G  R4
 L  R5
COLUMNS
    X1        COST         2.0         R1           1.0
    X1        R2           1.0         R4           1.0
    X1        R5           1.0
    X2        COST        -3.0         R1           1.0
    X2        R3           1.0         R5           1.0
    X3        COST         1.0         R2          -1.0
    X3        R5           1.0
    X4        COST         4.0         R1           1.0
    X4        R3           1.0         R4          -1.0
    X4        R5           1.0
    X5        COST         1.0         R2           1.0
    X5        R3           1.0         R4           2.0
RHS
    RHS       COST        -7.0
    RHS       R1          10.0         R2          -2.0
    RHS       R3           4.0         R4          -5.0
    RHS       R5          12.0
RANGES
    RNG       R4          13.0
BOUNDS
 FR BND       X1
 UP BND       X2           5.0
 FX BND       X3           2.0
 LO BND       X4          -3.0
 UP BND       X4           6.0
 LO BND       X5          -1.0
ENDATA
