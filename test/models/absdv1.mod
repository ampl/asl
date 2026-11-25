# abs in a defined variable, subtracted after a quotient and a constant:
# the solvers2 pfgh_read freed the adjoint cell of abs's argument twice,
# so the quotient's partial derivative overwrote the cell of the sum.
var x0 := 1.98;

var d0 = abs(1.64 - x0);

minimize Obj: (x0/(2 + x0^2) - 0.80) - d0;
