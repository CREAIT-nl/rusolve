\ G1 (e): LP for passLinearObjectives
minimize
 obj: x1 + x2 + x3 + x4
subject to
 r1: x1 + x2 + x3 >= 2
 r2: x2 + x3 + x4 >= 3
 r3: x1 - x4 <= 1
 r4: x1 + x2 + x3 + x4 <= 10
bounds
 0 <= x1 <= 5
 0 <= x2 <= 5
 0 <= x3 <= 5
 0 <= x4 <= 5
end
