\ G1 (b): -inf cost fixes x at its upper bound 3, violating c1
\ Infinite costs are written 1e400 (read as +inf): the HiGHS 1.16 LP
\ reader rejects the inf keyword as a coefficient.
minimize
 obj: - 1e400 x + y
subject to
 c1: x + y <= 2
bounds
 1 <= x <= 3
 0 <= y <= 1
end
