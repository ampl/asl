# shared defined variable whose body is another defined variable (f = d):
# its derivatives must reach d's own adjoint (fg_read)
var x := 1.3; var y := 2.1; var z := 0.7;
var d = x*y;
var f = d;
minimize o: f*f*z;
s.t. c2: f*z <= 40;
