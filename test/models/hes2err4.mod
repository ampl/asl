# As hes2err1, with the undefined second derivative (Q^1.852 at Q = 0)
# in the objective only: an error with the objective's weight, none when
# only the constraint has weight.
var Q >= 0, <= 0, := 0;
var y := 0.7;
var z := -0.3;

minimize o: (Q^1.852 + y)^2 + z^2;

s.t. c: y^2 + exp(z) + Q*y >= 1;
