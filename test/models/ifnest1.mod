# if in the else branch of an if whose then branch has derivatives
# (fg_read: the result cell must not be one the then branch uses)
var x := 1.3; var y := 2.1; var z := 0.7;
s.t. c1: (if y >= 0 then x^2 else (if 0.8 >= z then y else z)) <= 1e6;
s.t. c2: x + (if y >= 0 then x^3 else (if 0.8 >= z then y*z else 3)) <= 1e6;
