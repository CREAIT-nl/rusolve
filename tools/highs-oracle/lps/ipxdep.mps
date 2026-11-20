* rusolve highs-oracle synthetic LP for the IPX fixtures: linearly
* dependent equality rows (D3 = D1 + D2, D6 = D4 - D5) and duplicate
* columns (Z4 = Z5), so that the crash basis needs Basis::Repair and the
* starting basis reports dependent rows / columns.
NAME          IPXDEP
ROWS
 N  COST
 E  D1
 E  D2
 E  D3
 E  D4
 E  D5
 E  D6
 L  D7
COLUMNS
    Z1        COST         1.0         D1           1.0
    Z1        D3           1.0         D4           1.0
    Z1        D6           1.0         D7           1.0
    Z2        COST         2.0         D1           1.0
    Z2        D3           1.0         D5           1.0
    Z2        D6          -1.0
    Z3        COST         3.0         D1           1.0
    Z3        D2           1.0         D3           2.0
    Z3        D4           2.0         D6           2.0
    Z4        COST         1.0         D2           1.0
    Z4        D3           1.0         D5           1.0
    Z4        D6          -1.0         D7          -1.0
    Z5        COST         1.0         D2           1.0
    Z5        D3           1.0         D5           1.0
    Z5        D6          -1.0         D7          -1.0
    Z6        COST        -1.0         D2           1.0
    Z6        D3           1.0         D7           1.0
RHS
    RHS       D1           4.0         D2           6.0
    RHS       D3          10.0         D4           5.0
    RHS       D5           3.0         D6           2.0
    RHS       D7           8.0
BOUNDS
 UP BND       Z6           4.0
ENDATA
