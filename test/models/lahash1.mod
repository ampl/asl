# The linear-argument hash table grows while d1 is read.  The reader used
# to cut a hash chain after the rehash: the second use of x[1] (in d2) got
# its own linear argument, the objective's range held x[1] twice, and the
# Hessian kept by hvcomp (solvers2) had a wrong column for x[1].
# Any change in the number of linear arguments moves the rehash.
var x {1..6} >= 0.5, <= 3, := 1.3;

var d1 = x[1]*x[2] + 2*x[3] + x[5]^2 + 0.5*x[6] + x[3]*x[6]/5;
var d2 = log(x[1]);
var d3 = d1*d2;
var d4 = x[4]^2;
var lin = x[4] + x[5];

minimize o: d3 + d4;

s.t. c5: d1^2 <= 500;
s.t. c7: lin^2 + d4 <= 80;
