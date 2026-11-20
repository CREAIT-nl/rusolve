* LP writer edge cases: a quadratic objective (QUADOBJ, lower triangle),
* a variable that appears only in the Hessian, and a linear constant.
NAME          EDGEQP
ROWS
 N  obj
 G  c1
COLUMNS
    x         obj       1            c1        1
    y         obj       -2           c1        1
    z         obj       0
RHS
    rhs       obj       -3           c1        1
QUADOBJ
    x         x         2
    x         y         1.5
    y         y         4
    z         z         0.5
    x         z         -1
ENDATA
