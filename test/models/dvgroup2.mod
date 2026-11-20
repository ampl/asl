# As dvgroup1, with the defined variable split into 10 segments over four
# variables: the scratch range overruns in both libraries.
set I := 1..4;
var x {i in I} := 0.6 + 0.2*i;

var d = sum {i in I} x[i]^3
      + sum {i in I, j in I: i < j} x[i]*x[j];
var e = cos(x[1]*x[2]);

minimize o: d*x[1];

s.t. c1: d + e^2 <= 1000;
