# Q is fixed at 0 by its bounds (as a branch-and-bound node does), and
# Q^1.852 has no second derivative there (1 < p < 2): sphes after
# sphsetup, duthes and fullhes must report pow", as hvcomp does.
var Q >= 0, <= 0, := 0;
var y := 0.7;
var z := -0.3;

minimize o: (Q^1.852 + y)^2 + y*z + z^2 + Q*z;

s.t. c: Q^1.852 * y + exp(z) + y^2 >= 1;
