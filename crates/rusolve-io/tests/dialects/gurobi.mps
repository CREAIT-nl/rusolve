NAME 
OBJSENSE MAX
ROWS
 N  OBJ
 L  c1      
 E  c2      
 E  c3      
 G  c4      
 G  c5      
COLUMNS
    MARKER    'MARKER'                 'INTORG'
    x         OBJ       3
    x         c1        1
    x         c2        1
    x         c3        1
    MARKER    'MARKER'                 'INTEND'
    y         OBJ       2
    y         c1        1
    y         c2        -1
    y         c4        1
    y         c5        -1
    z         OBJ       -1
    z         c1        1
    z         c5        1
    MARKER    'MARKER'                 'INTORG'
    b         OBJ       4
    b         c3        1
    MARKER    'MARKER'                 'INTEND'
    w         OBJ       0.5
    v         OBJ       1
    v         c4        1
    Rgc2      c2        1
RHS
    RHS1      OBJ       -7
    RHS1      c1        10
    RHS1      c2        8
    RHS1      c3        3
    RHS1      c4        -5
    RHS1      c5        -4
BOUNDS
* Note: variable lower bounds are 0 unless specified otherwise
 UP BND1      x         6
 LO BND1      y         -2
 UP BND1      y         5
 FR BND1      z       
 BV BND1      b       
 FX BND1      w         2
 MI BND1      v       
 UP BND1      v         3
 UP BND1      Rgc2      6
ENDATA
