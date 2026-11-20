* SCIP STATISTICS
*   Problem name     : model
*   Variables        : 6 (1 binary, 1 integer, 0 implicit integer, 4 continuous)
*   Constraints      : 6
NAME          model
OBJSENSE
  MAX
ROWS
 N  Obj 
 L  c1 
 G  c2a 
 L  c2b 
 E  c3 
 G  c4 
 G  c5 
COLUMNS
    INTSTART  'MARKER'                            'INTORG'                           
    x         Obj                              3  c3                               1 
    x         c2a                              1  c1                               1 
    x         c2b                              1 
    b         Obj                              4  c3                               1 
    INTEND    'MARKER'                            'INTEND'                           
    y         c1                               1  c2a                             -1 
    y         c2b                             -1  Obj                              2 
    y         c4                               1  c5                              -1 
    z         Obj                             -1  c5                               1 
    z         c1                               1 
    w         Obj                            0.5 
    v         Obj                              1  c4                               1 
RHS
    RHS       c1                              10  c2a                              2 
    RHS       c2b                              8  c3                               3 
    RHS       c4                              -5  c5                              -4 
    RHS       Obj                             -7 
BOUNDS
 BV Bound     b                                  
 UP Bound     x                                6 
 FR Bound     z                                  
 LO Bound     y                               -2 
 UP Bound     y                                5 
 FX Bound     w                                2 
 MI Bound     v                                  
 UP Bound     v                                3 
ENDATA