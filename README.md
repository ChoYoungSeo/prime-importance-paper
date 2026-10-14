# PRIME: A New Measure for Variable Importance

[Paper](https://arxiv.org/abs/2610.17053) · [R package](https://github.com/ChoYoungSeo/prime-importance-r) · [Python package](https://github.com/ChoYoungSeo/prime-importance) · [Examples](#quick-start)

## Why PRIME?


PRIME measures which inputs explain a response and how much each contributes, even when observations are noisy and inputs are correlated. It estimates **proportional marginal effects (PME)**: under the paper's identifiability conditions, redundant inputs receive zero population importance. This addresses a limitation of Shapley effects, which can give substantial credit to correlated proxies.

<p align="center"><img src="docs/assets/proxy-allocation.png" alt="Shapley and PME importance as a correlated proxy becomes active" width="700"></p>

**A correlated proxy.** In **m(X) = X₁ + X₂ + βX₃**, the inputs are jointly Gaussian with standard-normal marginals; X₁ and X₃ have correlation ρ with |ρ| < 1, and X₂ is independent. At β = 0, the regression mean does not depend on X₃: its Shapley effect is ρ²/4, while its PME is zero.

PRIME estimates these signal contributions from observed data, assigning zero outside the selected set. The paper establishes consistency under its estimation and selection conditions; restricting allocation to selected inputs also reduces the cost of exact computation.

## Results


**Friedman.** Five active inputs and five correlated proxies. In the paper’s Table 2, PRIME has the lowest mean PME estimation error among PRIME, PRIME-Full, and PME-NA across all nine sample-size and correlation settings. The plot shows 30 paired datasets with n = 100 and ρ = 0.9; PME-NA denotes noise-adjusted PME-kNN with six neighbors, and FIRST-SHAP and mean absolute SHAP provide Shapley-based comparisons. Mean absolute SHAP is shown in response units.

<p align="center"><img src="docs/assets/friedman-current.png" alt="Friedman variable importance: PRIME, PME-NA, FIRST-SHAP and mean absolute SHAP" width="700"></p>

**Concrete strength.** Across 100 repeated samples of 824 observations, PRIME gives concentrated importance to the water-to-binder ratio and age. All four methods use the same observations; mean absolute SHAP summarizes a fitted boosted-tree model.

<p align="center"><img src="docs/assets/concrete-current.png" alt="Concrete variable importance: PRIME, PME-NA, FIRST-SHAP and mean absolute SHAP" width="700"></p>

[Experiment scripts](experiments/) · [Inputs and saved results](results/numerical_20261011/) · [Numerical protocols](docs/numerical-evidence.md)

## Install

Choose the R or Python implementation for your own data.

### R

Requires R 4.2 or newer and a C++ build toolchain.

```r
install.packages("remotes")
remotes::install_github("ChoYoungSeo/prime-importance-r", upgrade = "never")
```

[R package guide](https://github.com/ChoYoungSeo/prime-importance-r) · [R reference manual](https://github.com/ChoYoungSeo/prime-importance-r/blob/main/docs/prime.pdf)

### Python

Requires Python 3.10 or newer. The optional plotting dependency is included below.

```sh
python -m pip install "prime-importance[plot]"
```

[Python package guide](https://github.com/ChoYoungSeo/prime-importance) · [PyPI](https://pypi.org/project/prime-importance/)

## Quick start

### R

```r
library(prime)
set.seed(5)
X <- matrix(runif(300 * 8, -1, 1), nrow = 300)
colnames(X) <- paste0("X", seq_len(ncol(X)))
y <- X[, 1] - X[, 2] + X[, 5]^2 + rnorm(300, sd = 0.08)
fit <- prime(X, y)
fit$selected_names
fit$effects
```

### Python

```python
import numpy as np
import prime
rng = np.random.default_rng(5)
X = rng.uniform(-1, 1, size=(300, 8))
y = X[:, 0] - X[:, 1] + X[:, 4] ** 2 + rng.normal(0, 0.08, size=300)
fit = prime.fit(X, y)
print(fit.selected_names)
print(fit.effects)
```

Each example selects `X1`, `X2` and `X5`; the other coordinates receive zero importance. The package guides describe diagnostic plots and other options.

## Citation

Cho, Y. and Joseph, V. R. (2026). *Measuring Variable Importance Using
Proportional Inclusive Marginal Effects*. [arXiv:2610.17053](https://arxiv.org/abs/2610.17053).
[DOI:10.48550/arXiv.2610.17053](https://doi.org/10.48550/arXiv.2610.17053).

[Licence and attribution](docs/licensing.md)
