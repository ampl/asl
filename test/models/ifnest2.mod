# As ifnest1, with conditions that switch between the hvseq points: both
# branches of the outer if and of the inner if are evaluated.
var x := 1.3; var y := 2.1; var z := 0.7;
s.t. c1: (if y >= 2.4 then x^2 else (if 0.8 >= z then y else z)) <= 1e6;
s.t. c2: x + (if y >= 2.4 then x^3 else (if 0.8 >= z then y*z else 3))
	<= 1e6;
