# defined variables with a term without derivatives, used linearly and
# nonlinearly, next to a constraint term without derivatives
var x := 1.3; var y := 2.1; var z := 0.7;
var d = x*y + max(0.5, floor(x));
var e = d*z;
minimize o: e*e + d;
s.t. c: e + d*d + min(1, floor(z)) <= 30;
