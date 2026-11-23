\ G1 (b): +inf and large finite costs on a minimization; z is fixed at
\ its lower bound 0, w (cost 1e15) is infinite only with infinite_cost=1e15
\ Infinite costs are written 1e400 (read as +inf): the HiGHS 1.16 LP
\ reader rejects the inf keyword as a coefficient.
minimize
 obj: x + 2 y + 1e400 z + 1e15 w
subject to
 c1: x + y + z + w >= 2
 c2: x - y + 2 z <= 3
bounds
 0 <= x <= 4
 0 <= y <= 4
 0 <= z <= 5
 1 <= w <= 3
end
