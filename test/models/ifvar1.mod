# if with a plain-variable branch and a constant branch, times another
# variable (pfgh_read: the variable's cell must not hold the result)
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: x*z + (if y >= 1 then y else 3)*x;
s.t. c: (if z >= 0 then y else 3)*x + (if y >= 3 then x else 2)*z <= 30;
