* General integers with default, infinite and negative bounds (the CPLEX LP
* default for a general integer is [0, +inf)).
NAME          EDGEGEN
ROWS
 N  obj
 L  c1
COLUMNS
    MARKER                 'MARKER'                 'INTORG'
    gpl       obj       1              c1        1
    gli       obj       1              c1        1
    gmi       obj       1              c1        1
    gfr       obj       1              c1        1
    gbig      obj       1              c1        1
    gneg      obj       1              c1        1
    MARKER                 'MARKER'                 'INTEND'
RHS
    rhs       c1        100
BOUNDS
 PL bnd       gpl
 LI bnd       gli       -5
 MI bnd       gmi
 FR bnd       gfr
 UI bnd       gbig      1e30
 LI bnd       gneg      -7
 UI bnd       gneg      -2
ENDATA
