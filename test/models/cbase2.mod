# A constant base 0 in elements (y * 0^(1 + z^2), ...): the derivatives
# of 0^v w.r.t. v are 0 for v > 0 (the exponents stay positive at every
# test point). solvers2 pfgh reported pow" for them (as cbase1 for
# group functions).
param a := 0;
var y := 0.5;
var z := 0.2;

minimize o: y * a^(1 + z^2) + a^(1 + y^2) * z + (y - 1)^2 + z^2;

s.t. c: y * a^(2 + y^2 + z^2) + y^2 + z <= 3;
