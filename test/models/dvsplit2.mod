# Defined variables used linearly in the objective are split into terms.
# A term of d2 is the defined variable d1 itself; d3 uses d2 nonlinearly.
# solvers pfgh_read lost the derivatives of d2^2 through d1: for such a
# term it remapped an unrelated derivative cell.
var x {1..3} >= 0.5, <= 3, := 1;

var d1 = x[1]*x[2];
var d2 = x[3]^2 + d1;
var d3 = d2^2 + d1;

minimize o: d3 + d2;
