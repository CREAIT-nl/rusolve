\ G1 (i): nearly parallel rows
minimize
 obj: x + y + z
subject to
 c1: x + y >= 1
 c2: x + 1.0001 y >= 1.0001
 c3: y + z >= 0.5
bounds
 x free
 y free
end
