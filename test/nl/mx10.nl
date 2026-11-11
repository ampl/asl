g3 0 1 0	# problem mx10
 3 1 1 0 0	# vars, constraints, objectives, ranges, eqns
 1 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 3 3	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 3 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
x3
0 1
1 2
2 0.5
r
2 -100
C0
o0
o2
v1
o11
2
v0
v2
o5
v0
n2
O0 0
o54
3
o2
v0
o11
3
v2
v1
v0
o2
v2
o12
2
v0
v1
o5
v2
n3
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
