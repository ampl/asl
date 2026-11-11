var x := 1.3; var y := 2.2; var z := 0.6; var d = if z > 0 then y else x; minimize o: d*d; s.t. c: d^3 >= -100;
