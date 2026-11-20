# a term without derivatives (min of a constant and floor) next to other
# terms in the same variables: it must contribute nothing to the Hessian
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: x*y*z;
s.t. c: x*y + min(1, floor(z)) <= 30;
