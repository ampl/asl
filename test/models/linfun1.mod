# A linear defined variable d used nonlinearly by several functions and in
# another defined variable: solvers pfgh_read gives it a linear Hessian
# funnel (P.linhesfun > 0), which no other fixture builds.
var x {1..4} >= 0.5, <= 3, := 1;

var d = x[1] + 2*x[2] + 3;
var e = x[3]*x[4] + d;

minimize o: d^2*x[3] + e^2;

s.t. c1: exp(0.1*d)*x[4] <= 100;
s.t. c2: d*x[3]^2 + e*x[1] <= 50;
s.t. c3: sin(d)*cos(x[4]) + e^3 <= 50;
