# A defined variable with several nonlinear terms (split into segments),
# used nonlinearly in the objective and linearly in a constraint.
var x := 1.3; var y := 2.1; var z := 0.7;
var d = x^2*y + exp(z)*y;
minimize o: d*x;
s.t. c: d + z^2 <= 100;
