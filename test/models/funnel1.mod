# A defined variable shared by the constraints only, with a Hessian funnel
# (hvseq mode h: hvcomp for the objective alone, then for the Lagrangian).
var x := 1.3;
var y{i in 1..3} := 0.5 + 0.2*i;
var z := 0.8;

# costly, one argument, used in three ranges: a Hessian funnel
var v = sin(x)*exp(x)*x^3/(1 + x^2) + cos(x)^2*sqrt(1 + x^2)
      + x^2*log(1 + x^2) + exp(-x)*sin(2*x)/(2 + cos(x));

minimize O0: z^4 + z^2;

s.t. c{i in 1..3}: v*y[i]^2 + y[i]^3 <= 5;
