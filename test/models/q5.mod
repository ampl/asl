var x := 1.3; var y := 2.2; var z := 0.6; minimize o: sum{i in 1..3} (if x > i/2 then y else z)*x^i; s.t. c{i in 1..2}: max(x,y,z)*min(x*i,y) + (if z>i then x else y)*z >= -100;
