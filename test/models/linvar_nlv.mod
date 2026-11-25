# Written by AMPL; nl/linvar_nlv.nl is this model's NL file with the
# header's count of nonlinear variables in objectives lowered from 2 to 1,
# so that x[2] (column 1) is counted as linear although it enters d,
# which the objective uses nonlinearly, through d's linear part.  AMPL has
# written such headers (see README).  The readers must still use the value
# of x[2] when they evaluate d.
var x {1..2} >= 0.5, <= 3, := 1.3;

var d = 1.5*x[2] + x[1]^2;

minimize o: d*x[1];
