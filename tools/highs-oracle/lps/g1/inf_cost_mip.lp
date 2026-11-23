\ G1 (b): infinite cost on an integer column of a MIP (x fixed at
\ ceil(0) = 0); w (cost 1e15) is infinite only with infinite_cost=1e15
\ Infinite costs are written 1e400 (read as +inf): the HiGHS 1.16 LP
\ reader rejects the inf keyword as a coefficient.
maximize
 obj: - 1e400 x + 2 y + 3 z + 1e15 w
subject to
 c1: x + y + z <= 4.5
 c2: y - z >= -1
 c3: w + z <= 3
bounds
 0 <= x <= 3.5
 0 <= y <= 3
 0 <= z <= 2.5
 0 <= w <= 1
general
 x y
end
