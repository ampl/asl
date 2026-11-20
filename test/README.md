# ASL numerical regression tests

Small `.nl` fixtures are read through every reader of both libraries
(`solvers`: `fg`, `fgh`, `pfgh`; `solvers2`: `fg`, `pfgh`) and checked by
`derivcheck`:

- gradients and Jacobians against central differences of the function values;
- Hessians (`hvcomp`, and `sphes`/`fullhes` where available) against
  differences of the Lagrangian gradient and against each other;
- function values, gradients, Jacobians and evaluation errors of each reader
  against the first reader's, record for record;
- values and derivatives again after the perturbed evaluations (stale state).

Each fixture is evaluated at its start point and at a perturbed point, in two
call orders (values first, derivatives first). Both points are moved strictly
inside the variable bounds.

A fixture with a `nl/<fixture>.err` file is also evaluated exactly at its
start point (test point `0e`, `derivcheck` order `e`): values first, then,
after an evaluation at a nearby point, derivatives first.
The file lists the routines that must report an evaluation error there
(nonzero `*nerror`), one `routine index` line each (`objval`/`objgrd` with
the objective index, `conval`/`jacval` with 0); every other routine must
succeed, and its results are still checked.

`hvseq` reads each fixture with `pfgh_read` (both libraries), goes through 12
points near the start point with a single `hvcomp` call at each, and then
compares every product with central differences of the Lagrangian gradient.
The points are kept strictly inside the variable bounds, the direction is
scaled with the point, and the difference step is reduced where a bound is
close. A point with an evaluation error is skipped; the test fails (exit
status 4) when fewer than two points could be compared, and a non-finite
product or difference is a mismatch.
This finds Hessian-vector products computed with the Hessian of an earlier
point, which `derivcheck` cannot see: it calls `hvcomp` once per variable at
each of its two points.

`hvseq` also runs in other modes (third argument), for the solvers2
library; `m` and `h` also for the solvers library:
- `z` calls `hvcomp` with `ihd_limit = 0` (the Hessian is not kept,
  `hvp0comp`);
- `d` and `s` sum `hvcompd` or `hvcomps` over the objective and the weighted
  constraints (the sparse `hvcomps` result is scattered by its indices);
- `D` and `S` do the same, but without evaluating the functions first:
  `hvcompd` and `hvcomps` evaluate them;
- `m` calls `hvcompd` for the objective as the first call at each point,
  then evaluates the constraints and compares the `hvcomp` product;
- `o` evaluates the objective, calls `hvcomp` for it alone, and then adds
  `hvcompd` for each constraint without evaluating the constraints first;
- `h` evaluates the objective, calls `hvcomp` for it alone, then evaluates
  the constraints and compares the `hvcomp` product of the Lagrangian at the
  same point: the funnels must be recomputed when the shared defined
  variables of the constraints get their values (`funnel1`);
- `F` and `G` work as `D`, except at the first point: there the `hvcompd`
  calls are the first evaluation calls of the run (`F`: constraints first,
  then the objective; `G`: the objective first);
- `q` evaluates the objective, calls `hvcomp` for it alone, then `hvcomp`
  for the Lagrangian, which evaluates the constraints itself; `Q` is the
  mirror (constraints first); `r` and `R` are `q` and `Q` with
  `ihd_limit = 0`;
- `y` calls only `hvcomp` for the Lagrangian at each point, `Y` with
  `ihd_limit = 0`.

In modes `D`, `S`, `m`, `o`, `h`, `F` and `G` a first pass evaluates the
functions at every point (in `F` and `G`, after the first), only to skip the
points with evaluation errors. Modes `q`, `Q`, `r`, `R`, `y`, `Y` have no such
pass: every call at a point is a first call there, which is when shared
defined variables of functions not yet evaluated hold no values, and a
product must not visit them (`an1`, `an2`, `sg2`, `sg5`, `dvfirst1`..
`dvfirst3`, `funnel2`). They are not run for the `skip*` fixtures, which
have points with evaluation errors.

A fourth argument gives constraints weight 0 (bit i for constraint i): a few
tests use it to check that the defined variables of functions with weight 0
are not visited (`dvfirst1`, `dvfirst2`, `zw2`).

`hvcompd` and `hvcomps` leave x "known", so `hvseq` calls `xunknown()` after
them. In the solvers library the modes other than `m` and `h` still give
wrong products with defined variables and are not run.

`degcheck` (solvers2) checks on every fixture that a function which
`degree_ASL` calls quadratic is accepted by `mqpcheckv`: drivers that pass
quadratic terms to the solver (Knitro) stop otherwise.

## Running

```sh
cmake -S . -B build -DBUILD_ASL_TESTS=ON
cmake --build build
ctest --test-dir build
```

Single case: `build/bin/derivcheck-asl test/nl/mx1 pfgh 1 v` (see the usage
comment in [derivcheck.c](derivcheck.c)). The relative tolerance is 1e-4; the
finite-difference step is 1e-5 (1 + |x|).

## Fixtures

`nl/*.nl` are text NL files; `models/*.mod` are the AMPL models they were
written from with `option presolve 0, pl_linearize 0, auxfiles '';` and
`write g<name>;`. The two `lg*.nl` files have no model source;
`dvmerge2.nl` is `dvmerge1.nl` with `x^3 + d` in place of `x^3 - d`, the
power written first (AMPL writes `d` first).

| Fixtures | Subject |
|---|---|
| `pr34_pow15`, `pow_08_*`, `sq1`, `d1`..`d5`, `defvar_nl_02_*` | powers with constant exponents, defined variables |
| `asinh_*`, `th1`, `lg2`, `lg720` | inverse hyperbolic, `tanh`, logistic at large arguments |
| `i5`, `i11`, `i14`, `ifthen_nl_01_*`, `n1`, `n3`, `q1`, `q4`, `q5` | `if-then-else`, constant or plain-variable branches |
| `mx1`..`mx10`, `n4`, `n6`, `n7`, `q2` | `min`/`max`, plain-variable operands, nesting |
| `mxsel1`..`mxsel4` | `min`/`max` selecting an operand that is itself a `min`, `max` or `if` (clamps, also in a defined variable) |
| `mxvar1`, `mxvar2` | `min`/`max` over plain variables next to other terms in the same variables (sparse Hessian) |
| `mxcon1` | `min` selecting its constant operand (sparse Hessian) |
| `ifcst1` | `if` with a constant branch next to other uses of that constant |
| `ifvar1` | `if` with a plain-variable branch and a constant branch, times another variable (Hessian) |
| `ifnest1`, `ifnest2` | `if` in the `else` branch of an `if` whose `then` branch has derivatives; in `ifnest2` both conditions switch between the `hvseq` points |
| `mxnd1`, `mxnd2` | `min`/`max` without an operand that has derivatives, after other terms, also in defined variables used linearly |
| `dvcon1`, `dvif1`..`dvif3` | shared defined variables: next to a constraint without them, in `if` branches, with another defined variable as body; in `dvif3` the conditions switch between the `hvseq` points |
| `dvlin1`..`dvlin3` | a shared defined variable whose product has a factor with a linear argument (`sin(x+y)`, `(x+y)^2`) |
| `dvfirst1`..`dvfirst3` | a defined variable holding a `max`, shared by constraints only, of one constraint, or shared by objectives other than 0: a product for another function, as the first call of the run, must not visit it |
| `ifcell1` | an `if` whose branches differ in derivatives (`floor(y)` or `y`), times another term: the result needs its own cells |
| `dvchunk1` | many defined variables starting with a `min`/`max` without derivatives: a block of derivative records overflows inside one of them |
| `an1`, `an2` | a defined variable holding a `max`, shared by constraints only (objectives only), in the range of an objective (constraint) term |
| `sg2`, `sg5` | a split defined variable with a `max`, of objective 1 only (constraint 1 only) |
| `funnelc1`, `funnelo1` | a defined variable of one constraint (of the objective) with a Hessian funnel: the funnel must be recomputed after that function is evaluated |
| `funnel2` | a defined variable shared by the objectives only, with a Hessian funnel (mirror of `funnel1`) |
| `zw2` | `dvfirst1` with a linear constraint: weights on the linear constraint only |
| `skip1`, `skipall1` | evaluation errors at some (all) `hvseq` points: points are skipped; with none left, `hvseq` must fail |
| `obj2max2` | as `obj2max1`, with both objectives on the same variables (one range): the range is visited, its defined variable of objective 1 is not |
| `funnel1` | a defined variable shared by the constraints only, with a Hessian funnel |
| `obj2max1` | a second objective with a defined variable holding a `max`: products for objective 0 must not visit it |
| `elnd1`..`elnd5` | terms without derivatives (`min`/`max` of `floor`, `count`) next to other terms, also with defined variables (Hessian) |
| `dvmerge1`, `dvmerge2` | a term without defined variables, then a term in the same variable with a defined variable of that function only: the merged element must keep the defined variable (Hessian) |
| `dvgroup1`, `dvgroup2` | a defined variable split into 2 (10) segments, in the objective and in a constraint next to a group term `e^2`: the scratch of `pfgh_read` must not overrun into the group's variable markers (`dvgroup1` overran in the solvers2 library only) |
| `absdv1` | `abs` in a defined variable subtracted after a quotient and a constant: the adjoint cell of the argument of `abs` must be freed once (gradient, solvers2) |
| `dvsplit1` | a defined variable split into several terms, used nonlinearly and in a constraint (order of the backward pass) |
| `p1`, `p4`, `p5`, `s1`, `t7`, `t9` | piecewise-linear terms inside nonlinear expressions |
| `v1`, `v4`, `y1`, `y3`, `y4`, `q3` | defined variables holding `min`/`max`/`if`, single-use and shared |
| `evalerr_03_*`, `pr35_obj_sqrt` | evaluation errors in an objective next to constraints; `.err`: the objective gradient fails at the start point, the constraints do not |
| `owner_*` | defined-variable segments whose owner field contradicts the header: must be rejected by `fg_read`, `fgh_read` (solvers) and `fg_read`, `pfgh_read` (solvers2); the solvers `pfgh` reader ignores the field and accepts them |

## Limits

- An evaluation error skips the comparisons of the routine that reported it,
  and any error skips the Hessian checks. Errors are asserted only at `0e`
  points; at the other points they are compared between readers.
- Points on a kink of `abs`, `min`, `max` or a condition would give
  finite-difference mismatches; fixtures are chosen to avoid them.
- `numberof` in algebraic expressions is not covered; `count` only in `elnd4`.
- `pfg_read` of the solvers library is not covered.
- Sequences of points are covered for `hvcomp`, and in the solvers2 library
  also for `hvcompd` and `hvcomps`. Only objective 0 enters the Lagrangian
  (`hvseq` passes objective index 0): objectives 1 and later are not covered.
- Known open defects of the solvers2 library, not covered by fixtures:
  - `sphes` and `fullhes` with an `if` (even a single one, e.g.
    `minimize o: if y >= 0.3 then tanh(x) else sin(x*y);`) or with nested
    `min`/`max` give wrong values or crash;
  - `hvcomp` is wrong in some models with nested `min` in the objective
    next to a `min` in a constraint (right with `ihd_limit = 0`), and with
    an `if` in defined variables;
  - an `if` or a `max` in a defined variable that becomes a Hessian funnel:
    `hvcomp` gives wrong products or stops (`bad *o = 124 in hfg_fwd`,
    `bad *o = 165 ...`): the funnel code `hfg_fwd` lacks some ops.
    `hvcompd` and `hvcomps` do not use funnels and are right;
  - `floor` with defined variables: `pfgh_read` may stop with
    `bad e->op = 13 in ewalkg`, and some such models crash the
    Hessian-vector products;
  - `pfgh_read`: a defined variable split into several terms (used
    linearly in a function) can get wrong derivatives, e.g. when another
    defined variable that uses it enters a constraint nonlinearly.
