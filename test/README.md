# ASL numerical regression tests

Small `.nl` fixtures are read through every reader of both libraries
(`solvers`: `fg`, `fgh`, `pfgh`; `solvers2`: `fg`, `pfgh`) and checked by
`derivcheck`:

- gradients and Jacobians against central differences of the function values;
- Hessians (`hvcomp`, and `sphes`/`fullhes` where available) against
  differences of the Lagrangian gradient and against each other;
- function values, gradients and Jacobians of each reader against the first
  reader's;
- values and derivatives again after the perturbed evaluations (stale state).

Each fixture is evaluated at its start point and at a perturbed point, in two
call orders (values first, derivatives first).

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
`write g<name>;`. The two `lg*.nl` files have no model source.

| Fixtures | Subject |
|---|---|
| `pr34_pow15`, `pow_08_*`, `sq1`, `d1`..`d5`, `defvar_nl_02_*` | powers with constant exponents, defined variables |
| `asinh_*`, `th1`, `lg2`, `lg720` | inverse hyperbolic, `tanh`, logistic at large arguments |
| `i5`, `i11`, `i14`, `ifthen_nl_01_*`, `n1`, `n3`, `q1`, `q4`, `q5` | `if-then-else`, constant or plain-variable branches |
| `mx1`..`mx10`, `n4`, `n6`, `n7`, `q2` | `min`/`max`, plain-variable operands, nesting |
| `p1`, `p4`, `p5`, `s1`, `t7`, `t9` | piecewise-linear terms inside nonlinear expressions |
| `v1`, `v4`, `y1`, `y3`, `y4`, `q3` | defined variables holding `min`/`max`/`if`, single-use and shared |
| `evalerr_03_*`, `pr35_obj_sqrt` | evaluation errors in an objective next to constraints |
| `owner_*` | defined-variable segments whose owner field contradicts the header: must be rejected by `fg_read`, `fgh_read` (solvers) and `fg_read`, `pfgh_read` (solvers2); the solvers `pfgh` reader ignores the field and accepts them |

## Limits

- Evaluation errors at a test point skip the affected comparisons; error
  counts themselves are not checked.
- Points on a kink of `abs`, `min`, `max` or a condition would give
  finite-difference mismatches; fixtures are chosen to avoid them.
- `count` and `numberof` in algebraic expressions are not covered.
- `pfg_read` of the solvers library is not covered.
- Sequences of points are covered for `hvcomp` only.
