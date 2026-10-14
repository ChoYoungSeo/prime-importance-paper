# Numerical evidence

## Current four-panel comparisons

The current four-panel graphics compare PRIME, PME-NA, FIRST-SHAP and mean
absolute SHAP on 30 Friedman datasets and 100 Concrete subsamples. All methods
use the same observations. PRIME, PME-NA and FIRST-SHAP use standardized
predictors; SHAP uses the original predictor values and is reported in response
units.

`results/numerical_20261011/` contains the paired inputs, saved SHAP values,
cross-validation results, timing records and figure checks. The original
three-method estimates remain in `results/numerical_20261001/`.

```sh
make validate-release
make artifacts
```

Public validation checks scientific input hashes, saved allocation vectors,
normalization and ratios, the nine-condition numeric comparison below, and
retained approved research graphics. Paper prose, paper source, full paper PDFs
and typeset table source are `NOT_IN_PUBLIC_EXPORT`; their checks are not
reported as passed. Figure generation uses R with ggplot2, patchwork and jsonlite
and writes comparison PDF and PNG graphics only to
`build/current_numerical_20261011/`, using saved data without refitting models.

## Archived PME comparison

The archived experiments compare PRIME, PRIME-Full (local-linear PME using all
coordinates), and PME-NA (the original `sensitivity` 1.27.1 PME-kNN routine
with six neighbors and `noise=TRUE`). FIRST-SHAP appears in the allocation
figures as a Shapley-effect comparison. Each procedure receives the same
predictor columns centered by their sample means and divided by their sample
standard deviations; responses and row order are unchanged. Internal rescaling
is disabled. Native signed estimates are retained.

### Friedman function

There are 30 paired datasets in each of nine conditions: n = 100, 200, or 500;
p = 10; and dependence ρ = 0, 0.5, or 0.9. Five inputs enter the regression
mean and five are paired dependent proxies. PRIME and PRIME-Full target
unit-signal PME ψ. PME-NA targets γψ on the response-variance scale, where
γ = Var{m(X)}/Var(Y) = 0.9597204012134519.

The error is half-L1 distance divided by the population target total:

```math
\mathrm{error}_m=\frac{1}{2c_m}\sum_{j=1}^{10}
\left|\widehat\theta_{mj}-c_m\psi_j\right|,
\qquad c_m=1\text{ (PRIME, PRIME-Full)},\quad c_m=\gamma\text{ (PME-NA)}.
```

FIRST-SHAP has a different allocation target and is excluded from this table.
The nine-condition, three-method numeric comparison (Table 2 values) uses the
full-precision
[earlier comparison CSV](../results/numerical_20261001/friedman/earlier_four_method_table.csv):
`PRIME`, `Full LL` (PRIME-Full), and `PME-kNN` (PME-NA). Validation compares these
27 values to the saved estimates; it does not read a typeset table source. The
values displayed below are rounded to three decimal places.

| n | ρ | PRIME | PRIME-Full | PME-NA |
|---:|---:|---:|---:|---:|
| 100 | 0 | 0.101 | 0.140 | 0.644 |
| 100 | 0.5 | 0.107 | 0.143 | 0.404 |
| 100 | 0.9 | 0.103 | 0.134 | 0.748 |
| 200 | 0 | 0.067 | 0.104 | 0.428 |
| 200 | 0.5 | 0.062 | 0.094 | 0.398 |
| 200 | 0.9 | 0.083 | 0.109 | 0.773 |
| 500 | 0 | 0.047 | 0.080 | 0.741 |
| 500 | 0.5 | 0.047 | 0.072 | 0.294 |
| 500 | 0.9 | 0.045 | 0.058 | 0.843 |

The earlier three-panel Friedman plot shows n = 100 and ρ = 0.9. Its 0.019 s, 1.706 s, and 0.029 s
labels use the later matched serial timing cohort: 30 calls each for PRIME,
PME-kNN with each of two noise settings, and FIRST-SHAP (120 timed calls,
four warm-ups excluded). Standardization is timed inside each pipeline;
loading and validation are excluded. The displayed PME-NA panel uses TRUE.
Thirteen coordinates lie outside the detail window [−0.30, 1.05].

### Concrete strength

All procedures use the same 100 samples of 824 distinct rows, drawn without
replacement from 1,030 observations, and eight predictors. Each sample is
standardized independently. The 2,400 plotted coordinates are retained without
clipping, projection or fitted-sum normalization of the comparator. PME-NA
has 314 negative coordinates; 29 coordinates lie outside the detail window.
Full-data estimates are archived as method-specific stability references;
The earlier Concrete plot displays no reference diamonds.

## Reproduction and provenance

From the repository root, `make validate-release` checks the public scientific
inputs and saved results for the earlier three-method and current four-panel
comparisons. `make artifacts` rebuilds research graphics from saved estimates
in `build/current_numerical_20261011/`; it does not fit models or generate a
paper.

Public manifests identify the native exports, scientific source/input receipts,
full-precision numeric comparison CSV and retained approved graphic assets.
Historical locations in receipts are provenance only; the public validator and
figure builder use repository-relative paths. Paper prose, full paper PDFs,
TeX/Bib sources and typeset table source are excluded from this public export
and recorded as `NOT_IN_PUBLIC_EXPORT`, rather than successful validations.
The frozen raw input RDS bundles are included with unchanged hashes. Historical
receipts preserve their recorded release provenance. Public workflow manifests
record this export scope separately from the frozen scientific source, inputs
and saved vectors.

### Archived replay source

Native R replay and SHAP refit source remain available for advanced inspection
of the original seeds, common standardization, model settings and protocols.
These interfaces were not validated as execution workflows for this public
export. The archived R runner's preflight still needs maintenance for selection
between historical and current presentation overlays. The supported public
workflow is the saved-result validation and research-graphic generation above.

An archived smoke receipt records nine native fits for one Friedman n=100,
ρ=0.9 dataset and one Concrete824 sample, with maximum native-vector discrepancy
below 5 × 10⁻¹⁶. It includes PRIME-Full for Friedman and both PME-kNN noise
options. This historical receipt does not certify the current public replay
interface or timing measurements.

FIRST's official supplement `code.zip` contains no separate LICENSE file;
its associated article identifies CC BY-NC 4.0. The supplement source is not
redistributed under this repository's MIT license. Its original import helper
and recorded protocols remain available for advanced source inspection.

A saved SHAP case can be inspected without refitting:

```sh
python3 experiments/15_current_numerical_20261011/rerun_shap.py \
  --case-id friedman_n100_rho0.9_r001 --check-only
```

The saved data and recorded settings are retained for inspection. No full
native-method rerun or SHAP refit is part of the public validation receipt.
