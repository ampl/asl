# if with a constant branch next to other terms: the result of the if must
# not be stored in the constant's cell, which every use of the constant shares
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: y*z + (if x >= 1.6 then 1 else 2)*y + x*x;
s.t. c: (if x >= 1 then x else 1)*y + y*z <= 30;
