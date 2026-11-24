Minimize
 obj: - 0.05 w0 - 0.08 w1 - 0.11 w2 - 0.14 w3 + [ 0.4 w0 ^ 2 + 0.12 w0 * w1 + 0.04 w0 * w2 + 0 w0 * w3 + 0.9 w1 ^ 2 + 0.24 w1 * w2 + 0.2 w1 * w3 + 1.6 w2 ^ 2 + 0.6 w2 * w3 + 2.5 w3 ^ 2 ] / 2
Subject To
 budget: w0 + w1 + w2 + w3 = 1
 cardinality: y0 + y1 + y2 + y3 <= 2
 hi_0: w0 - y0 <= 0
 hi_1: w1 - y1 <= 0
 hi_2: w2 - y2 <= 0
 hi_3: w3 - y3 <= 0
 lo_0: w0 - 0.1 y0 >= 0
 lo_1: w1 - 0.1 y1 >= 0
 lo_2: w2 - 0.1 y2 >= 0
 lo_3: w3 - 0.1 y3 >= 0
Bounds
 0 <= w0 <= 1
 0 <= w1 <= 1
 0 <= w2 <= 1
 0 <= w3 <= 1
Binaries
 y0 y1 y2 y3
End
