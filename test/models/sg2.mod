# A split defined variable with a max, of objective 1 only (a product for
# objective 0 alone must not visit its segments).
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = x1^2*x2 + max(x2, x3)*x3;
minimize O0: x2*x3^2 + x1^2;
minimize O1: 2*v1 + x1^2;
