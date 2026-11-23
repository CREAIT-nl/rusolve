\ SCIP STATISTICS
\   Problem name     : model
\   Variables        : 6 (1 binary, 1 integer, 0 implicit integer, 4 continuous)
\   Constraints      : 6
Maximize
 Obj: +4 b +3 x -1 z +2 y +0.5 w +1 v +7
Subject to
 c1: +1 x +1 y +1 z <= +10
 c2a: +1 x -1 y >= +2
 c2b: +1 x -1 y <= +8
 c3: +1 x +1 b = +3
 c4: +1 y +1 v >= -5
 c5: +1 z -1 y >= -4
Bounds
 0 <= b <= 1
 0 <= x <= 6
 z free
 -2 <= y <= 5
 2 <= w <= 2
 -inf <= v <= 3
Binaries
 b
Generals
 x
End
