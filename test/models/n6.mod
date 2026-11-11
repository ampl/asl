var xx{1..2} := 1; var t = xx[2]; var x = xx[1]; var h = (5*t-2)^2/4 + (2*x-8)^2 + x*t; minimize Obj: xx[1]^4; s.t. C: h <= 100;
