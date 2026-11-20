# log and sqrt undefined at some hvseq points (x < 0): hvseq skips them
var x := 0.05; var y := 1;
minimize o: log(x)*y + y^2;
s.t. c: sqrt(x)*y <= 10;
