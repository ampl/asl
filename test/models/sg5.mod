# A split defined variable with a max, of constraint 1 only, next to an
# objective and a constraint 0 on the same variable pairs.
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = max(x1, x2)*x2 + max(x2, x3)*x3;
minimize O0: x2*x3^2 + x1^2*x2;
s.t. C0: x2^2*x3 + x1*x2 <= 10;
s.t. C1: v1 + x1^3 <= 10;
