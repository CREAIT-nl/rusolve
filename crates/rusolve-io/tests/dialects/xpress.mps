* Written from user model (presolvestate=0x680061)
NAME          t
OBJSENSE  MAXIMIZE
ROWS
 N  __OBJ___
 L  c1      
 L  c2      
 E  c3      
 G  c4      
 G  c5      
COLUMNS
    x         __OBJ___  3
    x         c1        1
    x         c2        1
    x         c3        1
    y         __OBJ___  2
    y         c1        1
    y         c2        -1
    y         c4        1
    y         c5        -1
    z         __OBJ___  -1
    z         c1        1
    z         c5        1
    b         __OBJ___  4
    b         c3        1
    w         __OBJ___  0.5
    v         __OBJ___  1
    v         c4        1
RHS
    RHS00001  __OBJ___  -7
    RHS00001  c1        10
    RHS00001  c2        8
    RHS00001  c3        3
    RHS00001  c4        -5
    RHS00001  c5        -4
BOUNDS
 UI BND00001  x         6
 UP BND00001  y         5
 LO BND00001  y         -2
 FR BND00001  z       
 BV BND00001  b       
 FX BND00001  w         2
 UP BND00001  v         3
 MI BND00001  v       
ENDATA
