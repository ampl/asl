# A random model on which the linear-argument hash table of solvers
# pfgh_read (and solvers2) lost entries after a rehash, as in lahash1:
# a duplicated linear argument gave wrong gradients and Hessian-vector
# products.  Which model triggers it depends on when the table grows, so
# almost any change hides the defect.
var x {1..7} >= 0.5, <= 3, := 1.3;

var d0 = x[2]*x[5] + x[7]*x[5]/5 + x[3]^2;
var d1 = 0.5*x[1] + log(x[4]) + log(x[7]) + log(x[7]) + x[3]*x[1];
var d2 = 2*x[1] + 1.5*x[5];
var d3 = d1*x[6] + log(x[3]);
var d4 = x[6]^2 + x[7]*x[4];
var d5 = log(x[4]) + d2*x[2] + d4*x[1] + 2*x[5];
var d6 = d3*x[3] + x[7]^2 + x[6]^2 + d0*x[7];
# d7 is not used: AMPL does not write it
var d7 = 1.5*x[6] + x[3]*x[3]/5 + x[4]*x[4]/5 + x[7]*x[5] + d0*x[7];

minimize o: (d6 + d5)^2/10 + log(x[7]);

s.t. c0: (d6 + d5 + x[3])^2 <= 1e4;
s.t. c1: (d0 + d2 + x[2])^2 <= 1e4;
