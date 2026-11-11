var x := 1; var y := 2; minimize o: x*(if y > 0 then x else y) + y*(if x > 5 then x else y) + y^3;
