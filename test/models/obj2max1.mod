# A second objective uses a defined variable holding a max: a product for
# objective 0 alone must not visit it (it is not evaluated at the point).
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = max(x3, x2);
minimize O0: x1^2;
minimize O1: 2*v1;
