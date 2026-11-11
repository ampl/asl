option presolve 0;
var x >= 0;
var y >= 0 := 1;
minimize f: sqrt(x);
subject to nl: x*y >= 0;
subject to lin {i in 1..4}: x + i*y >= 0;
write gpr35_obj_sqrt;
