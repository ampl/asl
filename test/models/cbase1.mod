# A constant base 0 (e.g. a variable fixed at 0 by presolve): 0^v is 0
# for v > 0, so its derivatives w.r.t. v are 0. Group form (an outer
# function of a linear term) and, for fg_read, the same powers in fg.
# (In an element, y * 0^z, solvers2 pfgh still reports pow" for the Hessian.)
param a := 0;
var y := 0.5;
var z := 0.2;

minimize o: a^(y + 2*z) + (y - 1)^2 + z^2;

s.t. c: a^(y + z) + y^2 + z <= 3;
