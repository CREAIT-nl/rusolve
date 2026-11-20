\ LP format for model comments
\ header comment line

\ blank lines and comment-only lines everywhere
Minimize \ objective keyword with trailing comment
 obj: 2 x + y \ objective terms
Subject To
 c1: x + y >= 3 \ a >= row
 c2: x <= 4 \ a <= row
Bounds
 y <= 10 \ an upper bound
End \ the end
