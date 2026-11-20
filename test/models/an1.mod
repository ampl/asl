# A defined variable holding a max, shared by two constraints only, in the
# range of an objective term: hvcomp for the objective alone, as the first
# call of the run, must not visit it.
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = max(x3, x2);
minimize O0: x2*x3^2 + x1^2;
s.t. C1: v1^2 <= 10;
s.t. C2: v1*x2 <= 10;
