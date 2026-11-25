# As hes2err1, with the undefined second derivative only in the
# constraint: with weight 0 for the constraint (sphcheck ... y0) the
# Hessian is the objective's (H[Q,Q] = 2), and no error may be reported.
var Q >= 0, <= 0, := 0;
var y := 0.7;
var z := -0.3;

minimize o: Q^2 + (y - 1)^2 + z^2 + Q*z;

s.t. c: Q^1.852 * y + z >= -5;
