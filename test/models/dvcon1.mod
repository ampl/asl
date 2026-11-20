# defined variable shared by the objective and a constraint; a later
# constraint without defined variables but with a max (fg_read: its chain
# ends must not get the derivatives of the objective's defined variable)
var x := 1.3; var y := 2.1; var z := 0.7;
var d = x*y;
minimize o: d;
s.t. c1: d*x <= 40;
s.t. c3: z + max(x, y) >= 0;
