# the same next to a defined variable used linearly
var x := 1.3; var y := 2.1; var z := 0.7;
var d = x*y;
minimize o: d*d;
s.t. c: d + min(1, floor(z)) <= 30;
