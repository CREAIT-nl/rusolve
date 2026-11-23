\ LP format for model implicit_coeffs
\ Implicit coefficients: bare x, negated -x, glued 3x, spaced 2 x.
Minimize
 obj: x - y + 3z + 2 w
Subject To
 c1: x + y + z + w >= 4
 c2: -x + 2y - 3z + w <= 1
End
