# A term without defined variables followed by a term in the same variable
# that uses a defined variable of this function only: pfgh_read merges them
# into one element, which must keep the defined variable (Hessian 6x - 2).
# test/nl/dvmerge2.nl is this model with x^3 + d, x^3 written first.
var x := 1.98;

var d = (1.64 - x)^2;

minimize o: x^3 - d;
