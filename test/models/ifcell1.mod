# An if whose branches differ in derivatives (floor(y) has none), times
# another term: the result of the if needs its own cells, not the
# one-cell temporary of the branch without derivatives.
var x := 1.3; var y := 1.6;
minimize o: (if x > 0 then floor(y) else y)*x^2;
