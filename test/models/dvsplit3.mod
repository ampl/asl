# dvsplit2 with the defined variables also in a constraint (shared, not
# split per function).
var x {1..3} >= 0.5, <= 3, := 1;

var d1 = x[1]*x[2];
var d2 = x[3]^2 + d1;
var d3 = d2^2 + d1;

minimize o: d3 + d2;

s.t. c: d2 + d3 <= 100;
