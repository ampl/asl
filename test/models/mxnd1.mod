# min/max without an operand that has derivatives, after other terms of
# the same expression (the derivatives of those terms must be kept)
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: (x*y + max(0.5, floor(x)))*z;
s.t. c: (y*z + min(5, floor(y)))*x <= 30;
