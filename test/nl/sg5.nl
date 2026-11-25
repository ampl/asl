g3 0 1 0	# problem sg5
 3 2 1 0 0	# vars, constraints, objectives, ranges, eqns
 2 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 3 3	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 6 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 1 0	# common exprs: b,c,o,c1,o1
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
C0
o0
o2
o5
v1
n2
v2
o2
v0
v1
V3 0 2
o0
o2
o12
2
v1
v0
v1
o2
o12
2
v2
v1
v2
C1
o0
v3
o5
v0
n3
O0 0
o0
o2
v1
o5
v2
n2
o2
o5
v0
n2
v1
k2
2
4
J0 3
0 0
1 0
2 0
J1 3
0 0
1 0
2 0
G0 3
0 0
1 0
2 0
