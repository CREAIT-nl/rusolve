\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  3 x + 2 y - z + 4 b + 0.5 w + v + 7 Constant
Subject To
 c1: x + y + z <= 10
 c2: x - y + Rgc2 = 8
 c3: x + b = 3
 c4: y + v >= -5
 c5: - y + z >= -4
Bounds
\ Note: variable lower bounds are 0 unless specified otherwise
 x <= 6
 -2 <= y <= 5
 z free
 w = 2
 -infinity <= v <= 3
 Rgc2 <= 6
 Constant = 1
Binaries
 b
Generals
 x
End
