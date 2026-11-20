# max selects its nested min operand (no constant anywhere)
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: max(min(x, y), z)*y*z;
s.t. c: x*y*z <= 30;
