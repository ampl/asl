g3 0 1 0	# problem zw2
 3 4 1 0 0	# vars, constraints, objectives, ranges, eqns
 3 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 2 2	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 10 2	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 1 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
x3
0 1.2
1 0.7
2 0.9
r
1 10
1 10
1 10
1 5
V3 0 0
o12
2
v1
v2
C0
o2
v0
v2
C1
o5
v3
n2
C2
o2
v3
v0
C3
n0
O0 0
o2
o5
v0
n2
v1
k2
3
6
J0 2
0 0
2 0
J1 2
1 0
2 0
J2 3
0 0
1 0
2 0
J3 3
0 1
1 1
2 1
G0 2
0 0
1 0
