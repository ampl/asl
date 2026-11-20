# A defined variable shared by the objective and a constraint: a product
# whose last factor has a linear argument (sin(x+y)).
var x := 1.3; var y := 0.6;
var d = sin(x)*cos(y)*sin(x+y);
minimize o: d;
s.t. c: d <= 100;
