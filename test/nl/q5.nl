g3 0 1 0	# problem q5
 3 2 1 0 0	# vars, constraints, objectives, ranges, eqns
 2 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 3 3	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 6 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
x3
0 1.3
1 2.2
2 0.6
r
2 -100
2 -100
C0
o0
o2
o12
3
v2
v1
v0
o11
2
v1
v0
o2
o35
o29
v2
n1
v0
v1
v2
C1
o0
o2
o12
3
v2
v1
v0
o11
2
v1
o2
n2
v0
o2
o35
o29
v2
n2
v0
v1
v2
O0 0
o54
3
o2
o35
o29
v0
n0.5
v1
v2
v0
o2
o35
o29
v0
n1
v1
v2
o5
v0
n2
o2
o35
o29
v0
n1.5
v1
v2
o5
v0
n3
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
