# a term without derivatives first in the objective
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: max(0.5, floor(x)) + x*y*z;
s.t. c1: y*z*z <= 40;
