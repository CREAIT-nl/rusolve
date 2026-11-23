g3 1 1 0	# problem unknown
 6 5 1 1 1 	# vars, constraints, objectives, ranges, eqns
 0 0 0 0 0 0	# nonlinear constrs, objs; ccons: lin, nonlin, nd, nzlb
 0 0	# network constraints: nonlinear, linear
 0 0 0 	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 1 1 0 0 0 	# discrete variables: binary, integer, nonlinear (b,c,o)
 11 6 	# nonzeros in Jacobian, obj. gradient
 2 1	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
C0	#c1
n0
C1	#c2
n0
C2	#c3
n0
C3	#c4
n0
C4	#c5
n0
O0 1	#o
n7
x0	# initial guess
r	#5 ranges (rhs's)
1 10	#c1
0 2 8	#c2
4 3	#c3
2 -5	#c4
2 -4	#c5
b	#6 bounds (on variables)
0 -2 5	#y
3	#z
4 2	#w
1 3	#v
0 0 1	#b
0 0 6	#x
k5	#intermediate Jacobian column lengths
4
6
6
7
8
J0 3	#c1
0 1
1 1
5 1
J1 2	#c2
0 -1
5 1
J2 2	#c3
4 1
5 1
J3 2	#c4
0 1
3 1
J4 2	#c5
0 -1
1 1
G0 6	#o
0 2
1 -1
2 0.5
3 1
4 4
5 3
