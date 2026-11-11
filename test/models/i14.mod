var b >= 0 <= 1 := .7; var x >=3, <=6 := 4; var y := 1; minimize Obj: y^2; s.t. c: sin(x) + (if b>0.5 then 2 else x^2) + y^3 <= 100;
