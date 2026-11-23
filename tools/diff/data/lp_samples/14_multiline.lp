\ LP format for model multiline
\ Objective and rows wrapped across lines.
Minimize
 obj: 2 x +
   3 y +
   1 z
Subject To
 c1: x +
   2 y +
   3 z >=
   12
 c2: 4 x - z <=
   8
End
