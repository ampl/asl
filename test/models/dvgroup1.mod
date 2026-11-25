# A defined variable split into two segments on different variables, used
# in the objective and in a constraint next to a group term (e^2): in a
# small model the scratch range of pfgh_read overran into the marker array
# of the groups, and the group lost variables (wrong Jacobian).
var x := 1.3;
var y := 2.1;
var z := 0.7;

var d = x^2*y + exp(z)*y;
var e = cos(z*x);

minimize o: d*x;

s.t. c1: d + e^2 <= 100;
