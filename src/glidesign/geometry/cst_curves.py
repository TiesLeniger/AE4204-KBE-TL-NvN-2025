import numpy as np
import scipy
from typing import Union

def _cst_matrix_shape(poly_order: int, x: np.ndarray) -> tuple[int, int]:
    return len(x), poly_order + 1

def _cst_matrix(poly_order: int, x: np.ndarray, n1: float = 0.5, n2: float = 1.0) -> np.ndarray:
    # x has to be a 1D numpy array
    x = np.array(x).flatten()
    # Class function
    cls = x ** n1 * (1 - x) ** n2

    cst_mat = np.zeros(_cst_matrix_shape(poly_order, x))
    # Shape functions
    for o in range(0, poly_order + 1):
        s = scipy.special.binom(poly_order, o) * x ** o * (1 - x) ** (poly_order - o)
        cst_mat[:, o] = cls * s

    return cst_mat

def fit_cst(x: np.ndarray, z: np.ndarray, poly_order: int) -> tuple[np.ndarray, float]:
    # remove trailing edge thickness
    z_te = z[-1]
    z_modified = z - x * z_te

    # create cst matrix and compute coefficients
    cst_mat = _cst_matrix(poly_order, x)
    coeff = np.linalg.lstsq(cst_mat, z_modified, rcond=None)[0]

    return coeff, z_te

def fit_cst_airfoil(x: Union[list, np.ndarray], z: Union[list, np.ndarray], poly_order: Union[int, list[int], tuple[int]]) -> tuple[list, list, float, float]:
    # Ensure x and z are 1D numpy arrays
    x = np.array(x).flatten()
    z = np.array(z).flatten()

    LE_idx = np.argmin(x)
    x_u = np.flipud(x[0:LE_idx + 1])
    z_u = np.flipud(z[0:LE_idx + 1])
    x_l = x[LE_idx:]
    z_l = z[LE_idx:]

    def normalize(array: np.ndarray) -> np.ndarray:
        return (array - np.min(array)) / (np.max(array) - np.min(array))
    
    x_u = normalize(x_u)
    x_l = normalize(x_l)

    if isinstance(poly_order, (list, tuple)):
        if len(poly_order) == 1:
            order_u = poly_order[0]
            order_l = poly_order[0]
        elif len(poly_order) == 2:
            order_u = poly_order[0]
            order_l = poly_order[1]
        else:
            raise ValueError("The number of specified polynomial orders can be at most 2; 1 for the upper side and one for the lower side of the airfoil")
    else:
        order_u = poly_order
        order_l = poly_order

    cst_u, z_te_u = fit_cst(x_u, z_u, order_u)
    cst_l, z_te_l = fit_cst(x_l, z_l, order_l)

    return cst_u.tolist(), cst_l.tolist(), z_te_u.tolist(), z_te_l.tolist() 