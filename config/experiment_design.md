# Locked Experiment Design

## 1. Coordinate recovery

The coordinate-recovery study evaluates the PRIME screening stage.
It uses Branin, Ackley-3, and Friedman-1 in `p=100`
dimensions; `n` is 100, 200, 500, 1000, or 2000. Predictor laws are
independent, spaced AR(1), and balanced correlated proxies. Each cell has 50
replications. The PRIME screen is compared with FIRST, LASSO, FOCI, and FORD.

FIRST, LASSO, FOCI, and FORD provide the competing selected sets. Current
PRIME-screen results use the BIC-gated ridge-start screen with querywise
forward and backward calibration.

## 2. Selection followed by PME

At `n=p=100`, the independent and balanced-proxy designs are used with
Branin, Ackley-3, and Friedman responses and 30 paired replications. Each
selector supplies a coordinate set. Every completed run is passed to the
same local-linear subset-value estimator and the same exact PME operator. This
isolates the downstream effect of coordinate selection.

The displayed downstream comparison uses the PRIME screen, FIRST, FOCI, and FORD. Exact
enumeration is used for their nonempty selected sets; nonpositive full-set
cases are counted separately and are not assigned an arbitrary statistical
loss. LASSO remains in the coordinate-recovery study but is not part of the
downstream allocation figure.

The coordinatewise detail uses Friedman under balanced proxies at `n=200`,
`p=100`, and 30 paired replications. The PRIME screen, FIRST, FOCI, and FORD are
followed by the same local-linear PME estimator.

## 3. PME estimation and correlated proxies

Two complementary comparisons are required.

1. A low-dimensional estimator calibration uses Branin (`p=8`) and Bratley
   (`p=10`) with paired inactive proxies. It compares PRIME,
   the two-nearest-neighbor subset-value estimator used by PME-kNN, an NCMD
   plug-in estimator, and
   local-linear diagnostic controls over `n` in 100, 250, and 500 and latent
   proxy correlation in 0, 0.5, and 0.9.
2. The complete high-dimensional proxy experiment retains the established
   manuscript setting: Branin and Bratley, `p=20`, `n` in 100, 200, and 500,
   latent proxy correlation 0.9, and 20 paired replications. Full-coordinate
   PME-kNN requires all `2^p` subset values and is therefore evaluated through
   the stated low-dimensional experiment. At `p=20` it would require
   1,048,575 nonempty coalition values, so the high-dimensional comparison
   uses complete screen-then-PME pipelines. In addition to FIRST and NNVS,
   LASSO, FOCI, and FORD are each followed by the PME-kNN allocation stage.

Accuracy is half-L1 distance to the normalized PME of the oracle regression
function. Inactive mass and correlated-proxy mass are reported separately.

## 4. Ambient-dimension scaling

The established benchmark grid is retained: Branin, Friedman-1, OTL, piston,
and borehole at `n=500` and `p` in 20, 50, 100, 250, and 500, with 30 paired
replications. Inputs are independent, so a mean-sufficient reduction leaves
the population PME unchanged. The comparison is among complete
screen-then-PME procedures using PRIME, FIRST, NNVS, LASSO, FOCI, and
FORD screens. Allocation error, exact-allocation completion, subset-value
calls, and elapsed time are retained at the task level. The main figure
reports allocation error for OTL, piston, and borehole only;
runtime and retained-dimension endpoints are supplementary diagnostics.

## 5. Observed-data applications and computational supplement

Concrete Compressive Strength, Airfoil Self-Noise, and Auto MPG each use 200
paired bootstrap samples for PME and five repeats of five-fold nested
cross-validation for prediction. The primary comparisons are screened versus
full-coordinate LL-PME allocation, coalition-value count, and selected versus
full-coordinate test RMSE. The Concrete analysis additionally uses an
engineering representation that replaces water by the water-to-binder ratio.
Data Twinning is evaluated as a representative-query approximation under the
current ridge-start screen. Backend parity and permuted-control experiments remain
provenance checks.

## Common rules

- All methods receive the same generated dataset within a replication.
- Complete-pipeline exact PME is evaluated for retained sets of at most 12
  coordinates, corresponding to at most 4,095 nonempty subset values. Larger
  selected sets are recorded as computationally incomplete and receive no
  allocation loss.
- Hyperparameters are fixed across response mechanisms and dependence laws.
- Raw allocations, projected games, selected sets, errors, and runtimes are
  checkpointed before any summary is produced.
- Oracle-coordinate calculations are diagnostics and are not labeled as
  competing procedures.
- Main-text claims are written only from validated main-profile outputs.
