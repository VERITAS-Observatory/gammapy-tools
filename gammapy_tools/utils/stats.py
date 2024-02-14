"""
Stats

A set of statistical functions
"""

import numpy as np


def get_cdf(
    arr: np.array, weights: np.ndarray = None, normalize: bool = True
) -> np.array:
    """Get the continious distribution function for a given array

    Parameters
    ----------
        arr (np.ndarray)                        - Array to get the CDF of.
        weights (np.ndarray)                    - Weights for each array element.
                                                  Defaults to None -> weights = 1
        norm (bool)                             - If the CDF should be normalized to 1.
                                                  Defaults to True
                                                  If false the cumulative distribution is returned

    Returns
    ----------
        cdf (np.ndarray)            - CDF of the `arr`

    """

    cdf = np.zeros(arr.shape)
    if weights is None:
        weights = np.ones(arr.shape)
    assert weights.shape == arr.shape

    cdf[0] = arr[0] * weights[0]
    for i in range(1, arr.shape[0]):
        cdf[i] = cdf[i - 1] + arr[i] * weights[i]

    if normalize:
        cdf /= cdf[-1]

    return cdf
