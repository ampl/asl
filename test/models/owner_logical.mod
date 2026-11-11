var b1 binary; var b2 binary; var x >=-19 <=32; var y >= -12 <= 5;
minimize Obj: b1 - x + 3*y;
var d = abs(x-y);
s.t. Con: b1 && b2 ==> x<=y;
s.t. Cone: d + x^2 <= 15;
