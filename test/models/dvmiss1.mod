# The NL file dvmiss1.nl is written from this model and then edited: the
# V segment of one defined variable is removed, as in a broken NL file
# that declares more defined variables than it provides.  Every reader
# must reject it (they left a null expression, and qp_read users such as
# mqpcheckv then crashed).
var x := 1.5;
var y := 0.5;
var v1 = x^2 + y;
var v2 = exp(x) * y;

minimize o: v1 + v2 * y;

s.t. c1: v1 * y <= 4;
s.t. c2: v2 + x * v1 >= 1;
