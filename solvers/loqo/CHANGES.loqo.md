# LOQO/AMPL Changelog

## 20201030
Relinked with and updated ASL2, which fixes some more minor problems affecting Hessian computations.

## 20201018
Relinked with an updated ASL2, fixing a bug affecting Hessian computations.

## 20201005
Relinked with an updated ASL2, fixing a possible problem with piecewise linear terms,

## 20200905
Relinked to use latest ASL version.

## 20191219
Relink to fix a rarely seen (and possibly harmless) bug.

## 20190908
Relink to ignore any LOGWAIT keywords in the ampl.lic or ampl.netlic file.

## 20190315
Fix bugs that occasionally affected sparsity computations and that affected "objrep" when several objectives can be adjusted.

## 20181221
Relink to fix a bug with piecewise-linear terms when "option pl_linearize 0;" is specified in AMPL.

## 20181221
Relink to fix possible trouble with complicated uses of more than one imported function.

## 20181210
Relink to fix possible trouble with "and" and "or" expressions in 64-bit binaries.

## 20181120
Relink to ignore HEARTBEAT lines in the ampl.lic file.

## 20180816
Relink to fix possible trouble with identifying quadratic objectives in large problems.

## 20180609
Relink to fix a fault with an example that used option presolve 0 (a bad idea).

## 20180525
Relink to compute tanh(x) and tanh'(x) for large |x| without complaint.

## 20180519
Relink to compute tanh(x) for large x without complaint.

## 20180402
Relink to fix a bug with nonlinear "if" expressions. Wrong gradients were possible.

## 20180314
Relink to fix a bug that gave error message "bad *o = ... in heswork".

## 20180302
Relink to fix error messages "Bad *o = 159 ..." or "... 127 ..." and to fix a bug (e.g., fault) with reading some large .nl files.

## 20180121
Relink with current ASL to fix a rarely seen bug.

## 20171129
Relink with current solver-interface library (to be safe).

## 20170801
Relink to fix a bug, introduced 20170511, with derivatives of abs().

## 20170619
Relink to fix several obscure bugs.

## 20170515
Relink to fix a glitch that caused an error message of the form "bad *o = ... in hfg_fwd".

## 20170511
Relink to fix a bug with defined variables shared by several constraints or objectives: under complicated conditions, it was possible for derivative evaluation errors to be ignored.
Change "version" keyword to a single-word phrase (no value assigned), as with various other solvers.

## 20160831
Relink to fix a bug in computing Hessians products when the same variable appears alone as the "then" or "else" part of two or more if-then-else expressions.

## 20160608
Relink to fix a bug with expressions of the form expr^num (with num a numeric constant) in "group partially-separable" contexts.

## 20160506
First version to be made available from AMPL Optimization.