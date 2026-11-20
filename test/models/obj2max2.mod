# As obj2max1, but both objectives have an element on the same variables
# (one range): the product for objective 0 alone must not visit v1.
var x1 := 1.2; var x2 := 0.7;
var v1 = max(x1^2, x2);

minimize O0: x1^2*x2;
minimize O1: 2*v1*x1;
