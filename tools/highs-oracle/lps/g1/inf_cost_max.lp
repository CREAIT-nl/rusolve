\ G1 (b): -inf and +inf costs on a maximization; x is fixed at its lower
\ bound 0, y at its upper bound 0; w (cost 1e15) is infinite only with
\ infinite_cost=1e15, fixed at its upper bound 3
\ Infinite costs are written 1e400 (read as +inf): the HiGHS 1.16 LP
\ reader rejects the inf keyword as a coefficient.
maximize
 obj: - 1e400 x + 1e400 y + z + 1e15 w
subject to
 c1: x + y + z + w <= 10
 c2: x - z >= -4
bounds
 0 <= x <= 6
 -2 <= y <= 0
 0 <= z <= 8
 0 <= w <= 3
end
