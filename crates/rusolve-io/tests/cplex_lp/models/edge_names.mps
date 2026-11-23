* LP writer edge cases: names CPLEX accepts as they are, names it cannot
* read (renamed), names that look like keywords or exponents.
NAME          EDGENAMES
ROWS
 N  obj
 L  e1
 G  end
 E  x(1,2)
 L  1row
 L  r:colon
 L  rr/slash
COLUMNS
    x(1,2)    obj       1            e1        1
    B0(998)   obj       2            end       1
    e1        obj       1            x(1,2)    1
    E.5       obj       1            e1        1
    eabc      obj       1            end       1
    e_        obj       1            1row      1
    a+b       obj       1            r:colon   1
    1x        obj       1            rr/slash  1
    x;y       obj       1            e1        1
    end       obj       1            end       1
    Free      obj       1            e1        1
    inf       obj       1            e1        1
    x,        obj       1            e1        1
    a_b       obj       1            e1        1
    a-b       obj       1            e1        1
    a!"#$%&'?@_`{}~.  obj  1         e1        1
    x/y       obj       1            e1        1
    x[1]      obj       1            e1        1
    x^2       obj       1            e1        1
    x|y       obj       1            e1        1
    Rge1      obj       1            e1        1
RHS
    rhs       e1        10           end       1
    rhs       x(1,2)    2            1row      5
    rhs       r:colon   5            rr/slash  5
RANGES
    rng       e1        4
BOUNDS
 UP bnd       a_b       3
ENDATA
