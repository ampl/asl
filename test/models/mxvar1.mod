# min/max over plain variables next to another body on the same variables
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: max(x, z)*y*y;
s.t. c: x*y*z <= 30;
