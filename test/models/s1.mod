var x := 1.3; var y := 2.2; var z := 0.6; minimize o: z^2; s.t. c: y * <<0,2; -1,1,3>> x + x*<<1; 2,5>> z >= -100;
