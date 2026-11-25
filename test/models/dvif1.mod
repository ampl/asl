# shared defined variables in the branch of an if whose other branch is
# a constant (fg_read: both branch chains must get their derivatives)
var x := 1.3; var y := 2.1; var z := 0.7; var w := 3;
var d = x*y;
var e = w*y;
var f = e^2;
minimize o: d + e*x;
s.t. c1: (if d >= 1 then d else 1)*z <= 40;
s.t. c2: f*(if x >= w then 1.96 else z) <= 1e6;
