var x := 1.3; var y := 2.2; var z := 0.6; minimize o: abs(x)*y + z*max(abs(x), y) + x*(if y > 0 then (if z > 0 then y else x) else z);
