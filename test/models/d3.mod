var xx{1..3} := 1; var h = (5*xx[2]-2)^2/4 + exp(xx[1])*xx[2]; minimize Obj: h; s.t. C: h <= 100; s.t. D: xx[3]^2 + h >= -3;
