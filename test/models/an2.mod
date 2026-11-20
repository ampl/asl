# A defined variable holding a max, shared by objectives 1 and 2 only, in
# the range of a constraint term: hvcomp for the constraints alone, as the
# first call of the run, must not visit it.
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = max(x3, x2);
minimize O0: x1^2*x2;
minimize O1: v1^2;
minimize O2: v1*x2;
s.t. C0: x2*x3^2 <= 10;
