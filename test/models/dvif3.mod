# As dvif1, with conditions that switch between the hvseq points: the
# constant branches are taken at some points.
var x := 1.3; var y := 2.1; var z := 0.7; var w := 0.5;
var d = x*y;
var e = w*y;
var f = e^2;
minimize o: d + e*x;
s.t. c1: (if d >= 2.7 then d else 1)*z <= 40;
s.t. c2: f*(if x >= 2*w then 1.96 else z) <= 1e6;
