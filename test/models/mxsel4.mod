# nested if selected by max; nested min selected as a later operand
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: max(if z >= 0 then x else y, z)*y*z;
s.t. c: max(z, min(x*y, 5))*y <= 30;
