# As hes2err1, with the power's argument an expression of a fixed and
# a free variable, both 0 at the start point.
var Q >= 0, <= 0, := 0;
var y >= 0, := 0;
var z := -0.3;

minimize o: (Q + y)^1.852 + y*z + z^2;

s.t. c: (Q + y)^1.852 + exp(z) + y^2 >= 0.5;
