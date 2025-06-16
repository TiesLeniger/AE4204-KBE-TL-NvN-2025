import numpy as np
import scipy
from typing import Union

def _cst_matrix_shape(poly_order: int, x: np.ndarray) -> tuple[int, int]:
    """
    Calculates the shape of the CST matrix based on the polynomial order and input array.

    Args:
        poly_order (int): The polynomial order for the CST.
        x (np.ndarray): The 1D numpy array representing the x-coordinates of the airfoil.

    Returns:
        tuple[int, int]: The shape of the CST matrix as a tuple (num_points, poly_order + 1).
    """
    return len(x), poly_order + 1

def _cst_matrix(poly_order: int, x: np.ndarray, n1: float = 0.5, n2: float = 1.0) -> np.ndarray:
    """
    Generates the CST matrix for a given set of x-coordinates using a specified polynomial order.

    Args:
        poly_order (int): The polynomial order for the CST.
        x (np.ndarray): The 1D numpy array representing the x-coordinates of the airfoil.
        n1 (float, optional): Exponent used in CST formulation for scaling the x-coordinates. Default is 0.5.
        n2 (float, optional): Exponent used in CST formulation for scaling the (1 - x) part. Default is 1.0.

    Returns:
        np.ndarray: The CST matrix, where each row corresponds to a set of shape functions.
    """
    # Ensure x is a flattened numpy array
    x = np.array(x).flatten()

    # Class function based on CST formulation
    cls = x ** n1 * (1 - x) ** n2

    # Initialize the CST matrix with appropriate shape
    cst_mat = np.zeros(_cst_matrix_shape(poly_order, x))

    # Loop over polynomial orders to compute shape functions
    for o in range(0, poly_order + 1):
        s = scipy.special.binom(poly_order, o) * x ** o * (1 - x) ** (poly_order - o)
        cst_mat[:, o] = cls * s  # Assign the shape function values to the matrix

    return cst_mat

def fit_cst(x: np.ndarray, z: np.ndarray, poly_order: int) -> tuple[np.ndarray, float]:
    """
    Fits a CST model to the given x and z coordinates using the specified polynomial order.

    Args:
        x (np.ndarray): The 1D numpy array of x-coordinates.
        z (np.ndarray): The 1D numpy array of z-coordinates (camber line).
        poly_order (int): The polynomial order for fitting.

    Returns:
        tuple[np.ndarray, float]: The fitted CST coefficients and the trailing edge thickness.
    """
    # Remove trailing edge thickness from z (last point of z array)
    z_te = z[-1]
    z_modified = z - x * z_te  # Modify z by subtracting the trailing edge component

    # Create CST matrix and compute coefficients using least squares fitting
    cst_mat = _cst_matrix(poly_order, x)
    coeff = np.linalg.lstsq(cst_mat, z_modified, rcond=None)[0]

    return coeff, z_te

def fit_cst_airfoil(x: Union[list, np.ndarray], z: Union[list, np.ndarray], poly_order: Union[int, list[int], tuple[int]]) -> tuple[list, list, float, float]:
    """
    Fits a CST model to the upper and lower surfaces of an airfoil.

    Args:
        x (Union[list, np.ndarray]): The 1D list or numpy array of x-coordinates.
        z (Union[list, np.ndarray]): The 1D list or numpy array of z-coordinates (camber line).
        poly_order (Union[int, list[int], tuple[int]]): The polynomial order(s) to use for the CST fit.

    Returns:
        tuple[list, list, float, float]: The CST coefficients for the upper and lower surfaces, 
                                          along with the trailing edge thickness for both surfaces.
    """
    # Ensure x and z are 1D numpy arrays
    x = np.array(x).flatten()
    z = np.array(z).flatten()

    # Find the index of the leading edge (LE)
    LE_idx = np.argmin(x)

    # Split the x and z arrays into upper and lower surfaces based on the leading edge
    x_u = np.flipud(x[0:LE_idx + 1])  # Flip x for the upper surface
    z_u = np.flipud(z[0:LE_idx + 1])  # Flip z for the upper surface
    x_l = x[LE_idx:]  # Lower surface x-coordinates
    z_l = z[LE_idx:]  # Lower surface z-coordinates

    def normalize(array: np.ndarray) -> np.ndarray:
        """
        Normalizes an array to the range [0, 1].

        Args:
            array (np.ndarray): The array to normalize.

        Returns:
            np.ndarray: The normalized array.
        """
        return (array - np.min(array)) / (np.max(array) - np.min(array))

    # Normalize the upper and lower surface x-coordinates
    x_u = normalize(x_u)
    x_l = normalize(x_l)

    # Handle cases where poly_order is a list or tuple
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

    # Fit CST for upper and lower surfaces
    cst_u, z_te_u = fit_cst(x_u, z_u, order_u)
    cst_l, z_te_l = fit_cst(x_l, z_l, order_l)

    return cst_u.tolist(), cst_l.tolist(), z_te_u.tolist(), z_te_l.tolist()