# As dvlin1, with the linear argument in the first factor.
var x := 1.3; var y := 0.6;
var d = sin(x+y)*sin(x)*cos(y);
minimize o: d;
s.t. c: d <= 100;
