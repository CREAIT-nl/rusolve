\ G1 (e): MIP for passLinearObjectives
maximize
 obj: 3 x + 2 y + z
subject to
 c1: 2 x + 3 y + z <= 12
 c2: x + y + 2 z <= 8
 c3: x - y >= -2
bounds
 0 <= x <= 5
 0 <= y <= 5
 0 <= z <= 4
general
 x y
end
