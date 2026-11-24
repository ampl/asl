g3 0 1 0	# problem lahash1
 6 2 1 0 0	# vars, constraints, objectives, ranges, eqns
 2 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 6 6 6	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 7 6	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 2 0 0 1 2	# common exprs: b,c,o,c1,o1
b
0 0.5 3
0 0.5 3
0 0.5 3
0 0.5 3
0 0.5 3
0 0.5 3
x6
0 1.3
1 1.3
2 1.3
3 1.3
4 1.3
5 1.3
r
1 500
1 80
V6 2 0
2 2
5 0.5
o54
3
o2
v0
v1
o5
v4
n2
o2
n0.2
o2
v2
v5
V7 0 0
o5
v3
n2
C0
o5
v6
n2
V8 2 2
3 1
4 1
n0
C1
o0
v7
o5
v8
n2
V9 0 3
o43
v0
V10 0 3
o2
v6
v9
O0 0
o0
v10
v7
k5
1
2
3
4
6
J0 5
0 0
1 0
2 0
4 0
5 0
J1 2
3 0
4 0
G0 6
0 0
1 0
2 0
3 0
4 0
5 0
