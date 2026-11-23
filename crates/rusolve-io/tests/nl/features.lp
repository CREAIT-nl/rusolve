\* Source Pyomo model name=unknown *\

min 
obj:
-3.5 ONE_VAR_CONSTANT
+2 i
-3 b
+1 f
-1 u
+0.25 w
+1 p
+4 q

s.t.

c_e_eq_:
+1 i
+1 b
+1 q
= 2

r_l_rng_:
+1 f
+1 u
-1 p
>= -1

r_u_rng_:
+1 f
+1 u
-1 p
<= 6

c_u_le_:
+1 f
+2 q
<= 3

c_l_ge_:
-1 i
+1 u
>= -7

c_l_cst_:
+1 w
+1 p
>= 2

c_u_dup_:
+2 i
-1 b
<= 9

r_l_box_:
+1 f
+1 u
+1 p
>= -20

r_u_box_:
+1 f
+1 u
+1 p
<= 20

bounds
   1 <= ONE_VAR_CONSTANT <= 1
   -5 <= i <= 5
   0 <= b <= 1
   -inf <= f <= +inf
   -inf <= u <= 4
   1.5 <= w <= 1.5
   0 <= p <= +inf
   -3 <= q <= 2
general
  i
binary
  b
end
