# The constraint's first derivative is undefined at Q = 0 (sqrt). With
# the objective given by its weight (nobj = -1, ow = {1}) and the
# constraint with weight 0, a Hessian call must not report the error.
var Q >= 0, <= 0, := 0;
var y := 0.7;

minimize o: (y - 1)^2 + y^4;

s.t. c: sqrt(Q) * y + y^2 <= 4;
