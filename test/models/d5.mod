var xx{1..3} := 1; var h = (5*xx[2]-2)^2/4 + exp(xx[1])*xx[2]; var g = sin(xx[3]*xx[1]) + h; minimize Obj: xx[2]^3; s.t. C: g <= 100; s.t. D: h + xx[1] >= -3;
