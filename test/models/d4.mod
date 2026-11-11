var xx{1..3} := 1; var h = (5*xx[2]-2)^2/4 + exp(xx[1])*xx[2]; var g = sin(xx[3]*xx[1]); minimize Obj: g; minimize O2: g; s.t. C: h <= 100; s.t. D: h*g >= -3; s.t. E: g <= 10;
