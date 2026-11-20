# count() in an algebraic constraint: a term without derivatives
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: x*y*z;
s.t. c: x*y + count {i in 1..2} (x <= 5*i) <= 30;
