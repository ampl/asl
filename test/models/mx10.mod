var x := 1; var y := 2; var z := 0.5; minimize o: x*min(x,y,z) + z*max(y,x) + z^3; s.t. c: y*min(z,x) + x^2 >= -100;
