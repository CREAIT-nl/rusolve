\ G1 (b): +inf cost on a free column (no lower bound to fix it at)
\ Infinite costs are written 1e400 (read as +inf): the HiGHS 1.16 LP
\ reader rejects the inf keyword as a coefficient.
minimize
 obj: 1e400 x + y
subject to
 c1: x + y >= 1
bounds
 x free
 0 <= y <= 3
end
