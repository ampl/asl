# A defined variable shared by the objective and the constraint that
# cannot be evaluated at the start point (sqrt of a negative number):
# sphsetup after the evaluation error must return, not end the run.
var x := -1; var y := 2;
var v = sqrt(x) + y^2;
minimize o: v + x^2*y;
s.t. c: v*x + y^3 <= 5;
