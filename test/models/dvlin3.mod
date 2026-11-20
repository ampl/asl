# As dvlin1, with a power of a linear argument; the constraint uses the
# defined variable nonlinearly.
var x := 1.3; var y := 0.6;
var d = sin(x)*cos(y)*(x+y)^2;
minimize o: 2*d;
s.t. c: d^2 <= 100;
