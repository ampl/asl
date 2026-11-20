# Many funnelkind-2 defined variables that start with a min/max without
# derivatives and use another defined variable: a derp chunk boundary
# (8190 derps) falls inside one of the bodies.
var x {1..3} := 0.4, <= 2;
var e = sin(x[1])*cos(x[2])*exp(x[3])*sin(x[2])*cos(x[1])*cos(x[3]);
var d {i in 1..120} = max(floor(x[1]) + 2, 3 - i)*e
   *prod {j in 1..40} (1 + 0.01*j*sin(x[1 + (i + j) mod 3]));
minimize o: sum {i in 1..120} d[i]^2/1e4;
s.t. c: sum {i in 1..120} d[i]*x[1 + i mod 3] <= 1e6;
