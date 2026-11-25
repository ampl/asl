# clamp held by a defined variable used by the objective and a constraint
var x := 1.3; var y := 2.1; var z := 0.7;
var d = min(max(x*y, 0.5), 5);
minimize o: (d - 2)^2 + z*z;
s.t. c: d*z <= 30;
