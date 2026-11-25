# A defined variable holding a max, used by one constraint only: a product
# for another function must not visit it before that constraint is
# evaluated (hvseq modes F and G).
var x1 := 1.2; var x2 := 0.7; var x3 := 0.9;
var v1 = max(x3, x2);
minimize O0: x1^2*x2;
s.t. C0: x1*x3 <= 10;
s.t. C1: v1^2*x1 <= 10;
