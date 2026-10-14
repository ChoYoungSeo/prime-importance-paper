# PRIME: A New Measure for Variable Importance

[Paper](https://arxiv.org/abs/2610.17053) · [Code and results](https://github.com/ChoYoungSeo/prime-importance-paper)

PRIME estimates proportional global variable importance from noisy data with
dependent predictors. It combines local-linear residual-variance estimation,
response-relevant variable selection, and proportional allocation on the
selected coordinates.

## Installation

```bash
python -m pip install prime-importance
```

To use `prime.plot(fit)`, install plotting support:

```bash
python -m pip install "prime-importance[plot]"
```

## Quick example

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

Use `prime.screen(X, y)` for selection diagnostics and `prime.plot(fit)` for
an importance plot.

## Citation

Cho, Y., and Joseph, V. R. (2026). *Measuring Variable Importance Using
Proportional Inclusive Marginal Effects*. [arXiv:2610.17053](https://arxiv.org/abs/2610.17053).
[DOI: 10.48550/arXiv.2610.17053](https://doi.org/10.48550/arXiv.2610.17053).

License: MIT. Authors: Youngseo Cho and V. Roshan Joseph.
