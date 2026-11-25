# A negative constant base with an integer exponent: (-2)^w is defined
# for integer w (here 2: 4); its derivative w.r.t. w is not (an error).
# solvers2 pfgh returned the value of that term for the whole function
# (4 instead of 8 for the objective, 4 instead of 8 for the constraint).
param b := -2;
var w := 2;

minimize o: b^w + w^2;

s.t. c: b^w * w <= 100;
