g3 0 1 0	# problem dvmiss1
 2 2 1 0 0	# vars, constraints, objectives, ranges, eqns
 2 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 2 2 2	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 4 2	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 2 1 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
x2
0 1.5
1 0.5
r
1 4
2 1
V2 0 0
o5
v0
n2
V3 0 0
o2
o44
v0
v1
V2 0 0
o5
v0
n2
C0
o2
v4
v1
C1
o0
v3
o2
v0
v4
O0 0
o0
v2
o2
v3
v1
k1
2
J0 2
0 0
1 0
J1 2
0 0
1 0
G0 2
0 0
1 1
