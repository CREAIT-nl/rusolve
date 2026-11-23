* LP writer edge cases: values CPLEX's decimal reader misreads when they
* are written the conventional way (17 digits, big and small magnitudes),
* -0, 1e30 and the largest double.
NAME          EDGENUM
ROWS
 N  obj
 L  r1
 G  r2
 E  r3
COLUMNS
    a         obj       0.12345678901234567   r1   0.40711116506242395
    a         r2        0.020160098014093186  r3   4.8523507186042204e-05
    b         obj       1e30                  r1   1e-30
    b         r2        3.0667406494635265e-09 r3  -0
    c         obj       1.7976931348623157e308 r1  2.2250738585072014e-308
    c         r2        1.2345678901234567e200 r3  9.87654321e-150
    d         obj       0.1                   r1   0.30000000000000004
    d         r2        123456789.12345679    r3   -2.5
RHS
    rhs       r1        0.33333333333333331   r2   -0
    rhs       r3        1e22
BOUNDS
 LO bnd       a         -0.7071067811865476
 UP bnd       a         3.141592653589793
 UP bnd       b         1e19
ENDATA
