# the same in defined variables used linearly (pfgh_read: a defined
# variable without operations that have derivatives)
var x := 1.3; var y := 2.1; var z := 0.7;
var d = x*y*z + max(0.5, floor(x));
var e = max(0.5, floor(y)) + x*y*z;
minimize o: d + e;
s.t. c: d + x*x <= 30;
