*SENSE:Maximize
NAME          t
ROWS
 N  OBJ
 L  c1
 G  c2a
 L  c2b
 E  c3
 G  c4
 G  c5
COLUMNS
    MARK      'MARKER'                 'INTORG'
    x         c1         1.000000000000e0
    x         c2a        1.000000000000e0
    x         c2b        1.000000000000e0
    x         c3         1.000000000000e0
    x         OBJ       -3.000000000000e0
    MARK      'MARKER'                 'INTEND'
    y         c1         1.000000000000e0
    y         c2a       -1.000000000000e0
    y         c2b       -1.000000000000e0
    y         c4         1.000000000000e0
    y         c5        -1.000000000000e0
    y         OBJ       -2.000000000000e0
    z         c1         1.000000000000e0
    z         c5         1.000000000000e0
    z         OBJ        1.000000000000e0
    MARK      'MARKER'                 'INTORG'
    b         c3         1.000000000000e0
    b         OBJ       -4.000000000000e0
    MARK      'MARKER'                 'INTEND'
    w         OBJ       -5.000000000000e-1
    v         c4         1.000000000000e0
    v         OBJ       -1.000000000000e0
RHS
    RHS       c1         1.000000000000e1
    RHS       c2a        2.000000000000e0
    RHS       c2b        8.000000000000e0
    RHS       c3         3.000000000000e0
    RHS       c4        -5.000000000000e0
    RHS       c5        -4.000000000000e0
BOUNDS
 UP BND       x          6.000000000000e0
 LO BND       y         -2.000000000000e0
 UP BND       y          5.000000000000e0
 FR BND       z       
 BV BND       b       
 FX BND       w          2.000000000000e0
 MI BND       v       
 UP BND       v          3.000000000000e0
ENDATA
