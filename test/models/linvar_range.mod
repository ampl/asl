# Written by AMPL; nl/linvar_range.nl is this model's NL file with the
# variable of the linear term in d's V segment (x[3], index 2) changed to 7,
# which is neither a variable nor a defined variable: every reader must
# reject the file (AMPL has written such indices, see README).
var x {1..3} >= 0.5, <= 3, := 1.3;

var d = 1.5*x[3] + x[1]*x[2];

minimize o: d*x[1];
