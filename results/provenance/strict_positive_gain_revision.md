# Strict positive-gain revision

## Decision

The empirical forward cutoff `0.05*sqrt(log(p)/n)` was removed.  The local
path now adds the coordinate with the largest estimated normalized
explained-variance gain while that gain is strictly positive.  The ridge seed,
neighborhood rule, local nugget, backward contrast, and PME operator are
unchanged.

The fixed-subset theory supplies consistency but not a nonasymptotic uniform
bound of the former cutoff's order.  Under the separated population-path
condition, the forward-coverage proof requires only that uniform subset-value
error vanish: if `4*epsilon_n < min(eta, delta)`, the selected empirical gain
is greater than `delta/2` and therefore positive.

## Paired pilot

The old and revised rules were compared on 180 paired `n=100`, `p=100`
selection datasets: Branin, Ackley-3, and Friedman under independent, spaced
AR(1), and balanced-proxy designs, with 20 replications per cell.  The same
data, ridge set, local path, and backward implementation were used.

| Rule | Mean candidate size | Active inclusion | Exact recovery | Mean false positives |
|---|---:|---:|---:|---:|
| Former empirical cutoff | 10.0 | 0.8833 | 0.8833 | 0.0000 |
| Strict positive gain | 10.9 | 0.8556 | 0.8556 | 0.0000 |

The final set changed in 5 of 180 runs.  All five were `n=100` Friedman runs
in which the strict path reached the 16-coordinate local-fit dimension limit;
the enlarged candidate set changed the sequential backward path and one or
more active coordinates were removed.  The strict path reached the dimension
limit in 35 of 180 runs.  These findings require a full PRIME-only result refresh
and prevent reuse of the former PRIME rows.

The row-level paired audit is `strict_positive_gain_pilot.csv`.
