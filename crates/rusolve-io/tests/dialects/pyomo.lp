\* Source Pyomo model name=unknown *\

max 
o:
+7 ONE_VAR_CONSTANT
+3 x
+2 y
-1 z
+4 b
+0.5 w
+1 v

s.t.

c_u_c1_:
+1 x
+1 y
+1 z
<= 10

r_l_c2_:
+1 x
-1 y
>= 2

r_u_c2_:
+1 x
-1 y
<= 8

c_e_c3_:
+1 x
+1 b
= 3

c_l_c4_:
+1 y
+1 v
>= -5

c_l_c5_:
-1 y
+1 z
>= -4

bounds
   1 <= ONE_VAR_CONSTANT <= 1
   0 <= x <= 6
   -2 <= y <= 5
   -inf <= z <= +inf
   0 <= b <= 1
   2 <= w <= 2
   -inf <= v <= 3
general
  x
binary
  b
end
