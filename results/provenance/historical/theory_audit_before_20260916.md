# PRIME theory audit

## Result-to-proof map

| Main result | Scope |
|---|---|
| Selected-set equivalence | Exact under a unique minimal mean-sufficient set |
| Noise-risk decomposition | Exact under conditional mean-zero response noise |
| Local-linear drift removal | Fixed projected dimension under Conditions L1--L3 |
| Fixed-subset consistency | Conditions L1--L4 with vanishing KNN radius |
| Forward coverage | Ridge-seed stability, gain separation, and uniform gain consistency |
| Active retention | Uniform deletion-gap consistency and a vanishing querywise buffer |
| PME stability | Uniform subset-risk consistency and a unique largest null subset |
| End-to-end consistency | Bounded selected dimension and uniform post-selection consistency |

Every named result in the main manuscript has a corresponding proof in the
Supplement. The local numerical nugget and Data-Twinning approximation are
recorded separately from the population theory.

## Query-calibrated backward reduction

## Candidate statistic

At a current set `S` of size `d`, all full and reduced fits use the same
neighborhood size `K`. For query `i` and deletion `j`, define

```text
Z_ij = local_risk_i(S without j; K) - local_risk_i(S; K),
Dhat_j = mean_i Z_ij,
scale_j = sd_i(Z_ij) / sqrt(q),
t_j = scale_j * sqrt(2 * log(max(2, d))),
Dtilde_j = max(Dhat_j - t_j, 0).
```

The procedure removes one coordinate with `Dtilde_j = 0`, recomputes all
neighborhoods and losses on the reduced set, and stops when every calibrated
loss is positive. The quantity `scale_j` is a querywise empirical fluctuation
scale. It is not used as an inferential standard error.

## The result required by PRIME

Let `A` be an active set and suppose the ridge-start forward path returns a
candidate set containing `A`. Assume empirical deletion losses are uniformly
consistent over the finitely many pairs visited by backward cleaning. Suppose
also that every active deletion has a population gap bounded below by
`Delta > 0` and that

```text
max_j t_j = o_p(1).
```

Then every active coordinate is retained with probability tending to one.
Indeed, an active empirical deletion loss converges to a value at least
`Delta`, while its threshold vanishes. Sequential recomputation preserves the
same argument at every visited superset of `A`.

The final set may contain inactive coordinates. This does not alter the
population target: PME exclusion assigns zero to those coordinates, and
superset invariance leaves the active-coordinate allocation unchanged.
Uniform consistency of the noise-adjusted subset values then gives consistency
of the reduced PME allocation.

## Why no overlap correction is needed for this theorem

Nearest-neighbor overlap matters if `scale_j` is interpreted as a standard
error or if the paper claims asymptotically exact removal of every null
coordinate. PRIME requires neither claim. The empirical scale obeys

```text
scale_j^2 <= mean_i(Z_ij^2) / (q - 1).
```

Thus bounded empirical second moments and `log(M_n)/q -> 0` make the maximum
threshold vanish without independence among the query contributions.
Bounded-overlap consistency controls `Dhat_j`. Together these facts are
sufficient for active retention and PME consistency on a sufficient superset.

Thus cross-fitting and network-HAC calibration would add computation to prove a
stronger exact-selection statement that is not needed for the allocation
method. Selection experiments should instead report active inclusion, false
positive rate, selected dimension, and final PME error separately.

## Adoption decision

Use ridge-start forward screening, sequential universal-threshold backward
cleaning, and the noise-adjusted PME game. State active retention and
sufficient-superset allocation consistency. Do not present the backward rule as
a hypothesis test and do not claim exact model-selection consistency.
