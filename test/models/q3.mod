var x := 1.3; var y := 2.2; var z := 0.6; var d = max(x,y); var e = if z > 0 then x else d; minimize o: d*z + e*e; s.t. c: e*y + d^2 >= -100;
