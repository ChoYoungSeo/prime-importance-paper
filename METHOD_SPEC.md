# PRIME: Locked Specification

## Population quantities

Let `P` denote all observed coordinates and define, for `S subset P`,

```
V(S) = Var{E(Y | X_S)}
L(S) = E{Var(Y | X_S)}
V(S) + L(S) = Var(Y).
```

For `Y = m(X) + epsilon` with `E(epsilon | X) = 0`,

```
L(S) = E(epsilon^2) + E[{m(X) - m_S(X_S)}^2].
```

The common noise term cancels from every deletion loss and selected-set
total-effect numerator. Finite-sample forward and backward comparisons are
regularized by paired querywise empirical scales.

Since `E(Y | X_S) = E{E(Y | X_P) | X_S}`, the total-effect set
function of the oracle regression mean `m_P(X) = E(Y | X)` can be
written on any mean-sufficient retained set `C` as

```
T_C(B) = {V(C) - V(C \ B)} / V(C)
       = {L(C \ B) - L(C)} / V(C),  B subset C.
```

PRIME applies the canonical nonnegative-game proportional-value extension to
`T_C` on the selected set and pads by zeros outside `C`.
If `A subset C` is the unique inclusion-minimal mean-sufficient set, then
`C \ A` is the greatest null coalition of `T_C`, containing every null coalition. The canonical proportional-value
extension assigns zero to `C \ A`, and its allocation on `A` is identical to the
PME computed directly from `T_A`. Thus sure screening is sufficient for the
population allocation; backward cleaning reduces the finite-sample and
computational burden.

## Estimator and FIRST-inspired selection

1. Standardize predictors.
2. Estimate `L(S)` by the average residual variance from centered
   local-linear ridge fits.
   Forward ranking uses the untruncated dimensionless score
   `U(S) = 1 - L(S)/Var(Y)`.
3. Use `K_n(d) = ceil(max(1.5 sqrt(n), 3d))`, capped at `floor(n/2)`.
4. Use a slope-only numerical nugget `lambda_loc = 1e-6` in local fits.
5. Rank coordinates by standardized global ridge coefficients. Search nested
   prefixes from zero through `ceil(1.5 log(max(2,p)))`, subject to local-fit
   feasibility, and select the initial prefix by

```
BIC(r) = n log{RSS(r) / n} + rank(r) log(n).
```

   The ridge and local numerical regularizers are both `1e-6` and have
   separate roles.
6. For nested sets `S subset T` evaluated with the same neighborhood size,
   define one paired querywise contrast and its empirical query scale:

```
H_i(S,T) = g_i(S) - g_i(T),
Gamma(S,T) = mean_i H_i(S,T),
omega(S,T) = sd_i{H_i(S,T)} / sqrt(n_Q),
Gamma_cal(S,T;m)
  = max{Gamma(S,T) - sqrt(2 log(max(2,m))) omega(S,T), 0}.
```

7. Starting from the BIC-selected ridge prefix `C`, rank each remaining
   coordinate by the untruncated raw gain `Gamma(C,C+j)`. Add the leader while
   `Gamma_cal(C,C+j;|C|+1)` is positive.
8. Sequentially remove one coordinate with
   `Gamma_cal(C-j,C;|C|) = 0`, recomputing every neighborhood and deletion loss
   after each removal.
9. Subtract the risk evaluated on the retained set as the estimated noise
   floor, take the positive part of each subset-risk difference, normalize the
   total-effect function, project it onto the closed convex set of monotone functions with
   fixed endpoints zero and one, and apply the extended proportional value.

The implementation accepts n >= 10 and p >= 1. It reports allocation as
unavailable when retained explained variance is at most
sqrt(machine epsilon) times the sample response variance. The default
proportional-value calculation treats normalized game values at most
sqrt(machine epsilon) as numerical zeros. The consistency statements refer
to the mathematically defined projection and canonical allocation before
these fixed software tolerances.

All full and reduced fits in a backward comparison use the same `K`.
Reduced neighborhoods are recomputed in the reduced coordinates.

Selection forms the paired difference before calibration.  Subtracting a
provisional noise floor from each intermediate risk and truncating the two
risks separately does not preserve their difference when either estimate
crosses that floor.  After a mean-sufficient set has been selected, the risk
evaluated on that set estimates the common noise floor used in the final PME
set function.

## Query approximation

All observations are query locations by default. Data Twinning can select
representative standardized `(X,y)` rows for the outer average;
all observations remain eligible as neighbors, and the same query set is
used in every comparison.

## Empirical scope

Coordinate recovery and every downstream allocation use the current
ridge-start screen. Benchmark functions, dependence designs, replications, and seeds are
held fixed for paired comparison. Previously computed Shapley allocations are
not reused as PME results.
