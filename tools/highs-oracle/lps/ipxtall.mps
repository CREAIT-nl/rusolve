* rusolve highs-oracle synthetic LP for the IPX fixtures: 8 rows > 2 * 3
* columns, so the default ("Lukas") strategy dualizes it; boxed (Y2) and
* free (Y3) columns.
NAME          IPXTALL
ROWS
 N  COST
 L  T1
 G  T2
 G  T3
 L  T4
 L  T5
 G  T6
 E  T7
 L  T8
COLUMNS
    Y1        COST        -1.0         T1           1.0
    Y1        T2           1.0         T4           1.0
    Y1        T5           2.0         T6          -1.0
    Y1        T7           1.0         T8           1.0
    Y2        COST        -2.0         T1           1.0
    Y2        T2          -1.0         T3           1.0
    Y2        T5           1.0         T7           2.0
    Y3        COST         1.0         T3           1.0
    Y3        T4           1.0         T5           1.0
    Y3        T6           1.0         T7           1.0
RHS
    RHS       T1           8.0         T2          -4.0
    RHS       T3           1.0         T4          10.0
    RHS       T5          20.0         T6          -6.0
    RHS       T7          12.0         T8           7.0
BOUNDS
 UP BND       Y2          10.0
 FR BND       Y3
ENDATA
