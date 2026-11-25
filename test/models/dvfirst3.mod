# A defined variable holding a max, shared by objectives 1 and 2: the first
# product of the run is for a constraint (hvseq mode F), and must not visit
# it before any objective is evaluated.
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = max(x3, x2);
minimize O0: x1^2*x2;
minimize O1: v1^2;
minimize O2: v1*x1;
s.t. C0: x1*x3 <= 10;
