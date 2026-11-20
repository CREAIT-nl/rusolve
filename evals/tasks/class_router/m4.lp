Minimize
 obj: ca + 3.1 b
Subject To
 demand: a + b >= 60
Bounds
 15 <= a <= 40
 0 <= b <= 50
 0 <= ca <= 400
Semi-continuous
 a
General Constraints
 pa: ca = PWL ( a ) : (0, 0) (15, 30) (25, 90) (40, 110)
End
