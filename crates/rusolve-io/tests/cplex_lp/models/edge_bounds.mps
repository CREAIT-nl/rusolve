* LP writer edge cases: every bound shape, ranged rows of each sense,
* an empty row, an unused column, zero costs, an objective constant, max.
NAME          EDGEBND
OBJSENSE
    MAX
ROWS
 N  profit
 L  lim
 G  need
 E  bal
 E  eqpos
 E  eqneg
 L  lrng
 G  grng
 G  emptyrow
 L  plain
COLUMNS
    free      profit    1.5          lim       1
    free      need      2
    fixed     lim       -1           bal       3
    neglo     need      1            eqpos     1
    uponly    eqneg     1            lrng      2
    minf      profit    -2           grng      1
    negup     profit    1            plain     1
    invert    profit    1            plain     1
    huge      profit    1            plain     1
    big20     profit    1            plain     1
    negzero   plain     1
    unused    profit    0
    zcost     lim       4            need      -1
    MARKER                 'MARKER'                 'INTORG'
    ibin      profit    3            lim       1
    igen      profit    1            need      1
    idef      profit    1            bal       1
    imark     profit    1            bal       1
    MARKER                 'MARKER'                 'INTEND'
    bvcol     profit    1            plain     2
    licol     profit    1            plain     1
RHS
    rhs       profit    -7.25
    rhs       lim       10           need      -3
    rhs       bal       4            eqpos     2
    rhs       eqneg     5            lrng      8
    rhs       grng      1            emptyrow  -3
    rhs       plain     100
RANGES
    rng       lrng      3            grng      -4
    rng       eqpos     2.5          eqneg     -1.5
BOUNDS
 FR bnd       free
 FX bnd       fixed     2.5
 LO bnd       neglo     -4
 UP bnd       uponly    9
 MI bnd       minf
 UP bnd       minf      6
 UP bnd       negup     -2
 LO bnd       invert    5
 UP bnd       invert    3
 UP bnd       huge      1e25
 LO bnd       big20     -1e20
 UP bnd       big20     1e20
 LO bnd       negzero   -0
 UP bnd       ibin      1
 LO bnd       igen      -3
 UP bnd       igen      12
 BV bnd       bvcol
 LI bnd       licol     2
 UI bnd       licol     7
ENDATA
