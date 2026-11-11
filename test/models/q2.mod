var x := 1.3; var y := 2.2; var z := 0.6; minimize o: max(x, min(y, z)) * z + min(max(x,y), z^2) * y;
