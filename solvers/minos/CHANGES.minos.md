# Summary of recent updates to MINOS for AMPL

## 20201030
Relinked with and updated ASL2, which fixes some more minor problems.

## 20201018
Relinked with an updated ASL2, fixing a bug affecting Hessian computations.

## 20190908
Relink to ignore any LOGWAIT keywords in the ampl.lic or ampl.netlic file.

## 20190318
Fix bugs with "objrep" when several objectives can be adjusted or when the adjusted objective must be scaled.

## 20181120
Relink to ignore HEARTBEAT lines in the ampl.lic file.

## 20180525
Relink to compute tanh(x) and tanh'(x) for large |x| without complaint.

## 20180519
Relink to compute tanh(x) for large x without complaint.

## 20180503
Relink to fix a bug (wrong gradients) with some uses of defined variables.

## 20171129
Relink to fix trouble with unconstrained problems whose objective is a defined variable.

## 20170803
Relink to fix possible trouble with objrep when the last constraint replaces the objective.

## 20170802
Relink to fix a bug, introduced 20170511, with (e.g.)

     var x; minimize o: if x < 3 then 3 else if x > 6 then 6 else x;

## 20170511
Relink to fix a bug with defined variables shared by several constraints or objectives: under complicated conditions, it was possible for derivative evaluation errors to be ignored.

## 20160329
Obscure bug fix: relink to fix a differentiation bug with the mod function.

## 20151026
Fix bugs with "objno" and "objrep" when the problem has several objectives.

## 20151014
Fix a fault with the following trivial example when solved with "option presolve 0;":

    var x; minimize o: x;
    s.t. c: x >= 42;

## 20150814
MacOSX binaries relinked to catch errors not reported via errno in evaluating some math functions.

## 20150630
Fix some possible trouble with a single-use license.

## 20150424
Fix a rarely seen licensing glitch.

## 20150217
Fix a possible fault with "objrep" on problems with nonlinear constraints and a linear objective; change objrep default to 3.

## 20141124
Relink for better handling of imported functions that report an inability compute derivatives.

## 20141029
Fix a possible fault on problems that only have constraints.

## 20141013
Relink macosx binaries so licenses can consider both hostname and local hostname.

## 20140828
Fix a glitch seen only on a bizarre MS Windows system that got eror message "Bad LOCAL_MGR = 0.0.0.0" with a floating license. Only the MS Windows bundles are updated. If you have not seen the "Bad LOCAL_MGR" message, you do not need this update. For MINOS, only the 32-bit MS Windows bundle is affected. With the updated minos.exe, invoking "minos -v" will show ASL(20140826).

## 20140313
Add keyword objrep: whether to replace

     minimize obj: v;

with

     minimize obj: f(x);

when variable v appears linearly in exactly one constraint of the form

     s.t. c: v >= f(x);

or

     s.t. c: v == f(x);

Possible objrep values:

## 0 = no

## 1 = yes for v >= f(x) (default)

## 2 = yes for v == f(x)

## 3 = yes in both cases.

## 20131023
Ignore case in MAC addresses during license checks (an issue rarely seen). When ending execution under a floating license, try to read a reply from the license manager to circumvent bug sometimes seen in MS Windows.

## 20131018
Relink to extend library renaming: if an imported-function library name has "_32" or "_64" before the final "." and fails to load (perhaps after changing "32" to "64" or vice versa, as appropriate), try omitting the "_32" or "_64".

## 20130320
Relink MS Windows versions to make automatic starting of ampl_lic work better on some versions of MS Windows (not XP). It is still recommended not to rely on automatically starting ampl_lic.

## 20120320
Adjust license-check in Linux versions for use with FreeBSD.

## 20120120
Relink to handle library names with '.' in a directory name but not in the basename (i.e., the name of the imported-function library).

## 20120117
Relink to simplify using a 64-bit minos with a 32-bit AMPL or vice versa when imported functions are involved (loaded from a *.dll file). For a 64-bit minos, if the library name involves '.' and the final '.' is preceded by "_32", change the "32" to "64". Otherwise, if the library fails to load and there is a '.' in the name, insert "_64" before the final '.'. (For 32-bit solvers, the rules are similar, with the roles of "32" and "64" reversed.)

## 20111229
Relink for use with single-user licenses.

## 20111107
Permit use of single-user licenses.

## 20111003
When processing ampl.lic, ignore new keywords for ampl.netlic.

## 20110527
Relink to permit a quoted "hostname" for MGR_IP in the ampl.lic file for a floating license.

## 20110426
Tweak license checker to correct a rare problem on MS Windows systems.

## 20110117
Mention "minos" in the "No license for this machine" message.