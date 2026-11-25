# As funnel1, but the defined variable is used by one constraint only
# (a c1 defined variable): its Hessian funnel is computed in hvcomp
# for the objective alone, before the constraint is evaluated.
var x := 1.3;
var y{i in 1..3} := 0.5 + 0.2*i;
var z := 0.8;

# costly, one argument, used in three ranges of one constraint
var v = sin(x)*exp(x)*x^3/(1 + x^2) + cos(x)^2*sqrt(1 + x^2)
      + x^2*log(1 + x^2) + exp(-x)*sin(2*x)/(2 + cos(x));

minimize O0: z^4 + z^2 + x^2;

s.t. c: sum{i in 1..3} (v*y[i]^2 + y[i]^3) <= 15;
