g3 1 1 0	# problem features_mp_text
 7 7 1 7 0	# vars, algcons, objs, ranges, eqns
 0 0	# nonlinear cons, objs
 0 0	# network cons: nonlinear, linear
 0 0 0	# nonlinear vars in cons, objs, both
 0 0 0 1	# linear network vars; funcs; arith, flags
 1 1 0 0 0	# discrete vars: binary, integer, nonlinear (b,c,o)
 17 7	# nonzeros in Jacobian, gradients
 3 1	# max name lengths: cons/objs, vars
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b	#7 bounds (on variables)
3
1 4
4 1.5
2 0
0 -3 2
0 0 1
0 -5 5
r	#7 ranges (rhs's)
4 2
0 -1 6
1 3
2 -7
2 2
1 9
0 -20 20
C0	#eq
n0
C1	#rng
n0
C2	#le
n0
C3	#ge
n0
C4	#cst
n0
C5	#dup
n0
C6	#box
n0
O0 0	#obj
n-3.5
k6	#intermediate Jacobian column lengths (cumulative)
3
6
7
10
12
14
J0 3
4 1
5 1
6 1
J1 3
0 1
1 1
3 -1
J2 2
0 1
4 2
J3 2
1 1
6 -1
J4 2
2 1
3 1
J5 2
5 -1
6 2
J6 3
0 1
1 1
3 1
G0 7
0 1
1 -1
2 0.25
3 1
4 4
5 -3
6 2
