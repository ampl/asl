# A defined variable shared by the objectives only, with a Hessian funnel
# (mirror of funnel1: hvcomp for the constraints alone, then for the
# Lagrangian after the objective is evaluated).
var x := 1.3;
var y{i in 1..3} := 0.5 + 0.2*i;
var z := 0.8;

var v = sin(x)*exp(x)*x^3/(1 + x^2) + cos(x)^2*sqrt(1 + x^2)
      + x^2*log(1 + x^2) + exp(-x)*sin(2*x)/(2 + cos(x));

minimize O0: sum{i in 1..3} (v*y[i]^2 + y[i]^3);
minimize O1: v*z;

s.t. c: z^4 + z^2 <= 5;
