Minimize
 obj: [ 0.04 w0 ^ 2 + 0.012 w0 * w1 + 0.004 w0 * w2 + 0 w0 * w3 + 0.09 w1 ^ 2 + 0.024 w1 * w2 + 0.02 w1 * w3 + 0.16 w2 ^ 2 + 0.06 w2 * w3 + 0.25 w3 ^ 2 ] / 2
Subject To
 budget: w0 + w1 + w2 + w3 = 1
 ret: 0.05 w0 + 0.08 w1 + 0.11 w2 + 0.14 w3 >= 0.08
Bounds
 0 <= w0 <= 1
 0 <= w1 <= 1
 0 <= w2 <= 1
 0 <= w3 <= 1
End
