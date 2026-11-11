var x := 1.3; var y := 2.2; var z := 0.6; minimize o: z*(if x > 0 then y else x) + (if y > 5 then z else x)^2 + y*(if x < 0 then x^2 else z);
