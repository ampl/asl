# A defined variable holding a max, shared by two constraints only: the
# first product of the run is for the objective (hvseq mode G), and must
# not visit it before any constraint is evaluated.
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = max(x3, x2);
minimize O0: x1^2*x2;
s.t. C0: x1*x3 <= 10;
s.t. C1: v1^2 <= 10;
s.t. C2: v1*x1 <= 10;
