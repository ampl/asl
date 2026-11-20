# A single if whose branches both have derivatives: sphes and fullhes set
# up their structures after the functions were evaluated at this x, and
# must evaluate them again (the setup used the cells of the taken branch).
var x := 1.3; var y := 0.6;
minimize o: if y >= 0.3 then tanh(x) else sin(x*y);
