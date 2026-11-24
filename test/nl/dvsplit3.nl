g3 0 1 0	# problem dvsplit3
 3 1 1 0 0	# vars, constraints, objectives, ranges, eqns
 1 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 3 3	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 3 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 3 0 0 0 0	# common exprs: b,c,o,c1,o1
b
0 0.5 3
0 0.5 3
0 0.5 3
x3
0 1
1 1
2 1
r
1 100
V3 0 0
o2
v0
v1
V4 0 0
o0
v3
o5
v2
n2
V5 0 0
o0
v3
o5
v4
n2
C0
o0
v5
v4
O0 0
o0
v5
v4
k2
1
2
J0 3
0 0
1 0
2 0
G0 3
0 0
1 0
2 0
