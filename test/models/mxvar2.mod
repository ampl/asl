var x := 1.3; var y := 2.1; var z := 0.7;
minimize o: max(x, y)*z + min(y, z)*x*x;
s.t. c: max(x, z)*y*y + x <= 30;
