# clamp min(max(e, lo), hi) selecting the nested max, after and before
# another term of the same body
var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: (min(max(x*y, 0.5), 5) - 2)^2 + z*z;
s.t. c: z*z + (min(max(x*z, 0.5), 5) - 2)^2 <= 30;
