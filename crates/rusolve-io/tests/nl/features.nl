g3 1 1 0	# problem unknown
 7 7 1 2 1 	# vars, constraints, objectives, ranges, eqns
 0 0 0 0 0 0	# nonlinear constrs, objs; ccons: lin, nonlin, nd, nzlb
 0 0	# network constraints: nonlinear, linear
 0 0 0 	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 1 1 0 0 0 	# discrete variables: binary, integer, nonlinear (b,c,o)
 17 7 	# nonzeros in Jacobian, obj. gradient
 3 1	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
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
x0	# initial guess
r	#7 ranges (rhs's)
4 2	#eq
0 -1 6	#rng
1 3	#le
2 -7	#ge
2 2	#cst
1 9	#dup
0 -20 20	#box
b	#7 bounds (on variables)
3	#f
1 4	#u
4 1.5	#w
2 0	#p
0 -3 2	#q
0 0 1	#b
0 -5 5	#i
k6	#intermediate Jacobian column lengths
3
6
7
10
12
14
J0 3	#eq
4 1
5 1
6 1
J1 3	#rng
0 1
1 1
3 -1
J2 2	#le
0 1
4 2
J3 2	#ge
1 1
6 -1
J4 2	#cst
2 1
3 1
J5 2	#dup
5 -1
6 2
J6 3	#box
0 1
1 1
3 1
G0 7	#obj
0 1
1 -1
2 0.25
3 1
4 4
5 -3
6 2
