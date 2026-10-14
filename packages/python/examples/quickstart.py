import matplotlib.pyplot as plt
import numpy as np

import prime


def make_data(n: int = 300) -> tuple[np.ndarray, np.ndarray]:
    index = np.arange(1, n + 1, dtype=float)
    frequencies = np.sqrt(
        np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29], dtype=float)
    )
    base = 2 * np.mod(np.outer(index, frequencies), 1) - 1
    rho = 0.9
    proxy_scale = np.sqrt(1 - rho**2)
    X = np.column_stack(
        [
            base[:, 0],
            base[:, 1],
            rho * base[:, 0] + proxy_scale * base[:, 2],
            rho * base[:, 1] + proxy_scale * base[:, 3],
            base[:, 4],
            base[:, 5],
            base[:, 6],
            base[:, 7],
        ]
    )
    y = X[:, 0] - X[:, 1] + X[:, 4] ** 2 + 0.08 * np.sin(index * np.sqrt(31))
    return X, y


if __name__ == "__main__":
    predictors, response = make_data()
    result = prime.fit(predictors, response)
    print(result)
    prime.plot(result)
    plt.show()
