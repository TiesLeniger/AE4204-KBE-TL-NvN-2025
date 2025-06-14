# Python native imports

# Python third party imports
import numpy as np

def generate_naca4(num_points: int, m: int, p: int, xx: int, spacing: str = 'cosine', closed_TE: bool = True) -> np.ndarray:
    """
    Generates a NACA 4-digit airfoil (e.g., NACA 2412) based on analytical expressions.

    Parameters:
        num_points (int): Number of points per surface (upper/lower).
        m (int): Maximum camber in percentage of chord (e.g., 2 for NACA 2412).
        p (int): Location of maximum camber in tenths of chord (e.g., 4 for NACA 2412).
        xx (int): Maximum thickness in percentage of chord (e.g., 12 for NACA 2412).
        spacing (str): Point spacing, either 'cosine' or 'constant'. Default is 'cosine'.
        closed_TE (bool): If True, TE is closed. If False, small open TE. Default is True.

    Returns:
        np.ndarray: Array of shape (2*num_points, 2), with concatenated upper and lower coordinates.
    """
    m /= 100
    p /= 10
    xx /= 100

    if spacing == 'cosine':
        beta = np.linspace(0.0, np.pi, num = num_points)
        x_c = (1-np.cos(beta))/2
    elif spacing == 'constant':
        x_c = np.linspace(0, 1, num = num_points)                               # Array of x coordinates, normalised w.r.t chord
    else:
        raise ValueError(f"spacing parameter must be 'cosine' or 'constant', got {spacing}")
    
    mask1 = (x_c >= 0.0) & (x_c < p)
    mask2 = (x_c >= p) & (x_c <= 1.0)

    y_camber = np.zeros_like(x_c)                                               # Array of y coordinates of the camberline, normalised w.r.t chord. Initialised as zeros
    dy_c_dx = np.zeros_like(x_c)
    if (m > 0 and p > 0):
        y_camber[mask1] = (m/p**2)*(2*p*x_c[mask1] - x_c[mask1]**2)                 # y_camber definition before position of maximum camber
        y_camber[mask2] = (m/(1-p)**2)*(1 - 2*p + 2*p*x_c[mask2] - x_c[mask2]**2)   # y_camber definition after position of maximum camber
        dy_c_dx[mask1] = (2*m/p**2)*(p - x_c[mask1])                                # y_camber definition before position of maximum camber
        dy_c_dx[mask2] = (2*m/(1-p)**2)*(p - x_c[mask2])                            # y_camber definition after position of maximum camber

    a4 = -0.1036 if closed_TE else -0.1015                                      # Coefficient that determines open or closed TE
    thickness_distribution = (xx/0.2)*(
        0.29698*np.sqrt(x_c) 
        - 0.126*x_c 
        - 0.3516*x_c**2 
        + 0.2843*x_c**3 
        + a4*x_c**4)

    theta = np.arctan(dy_c_dx)

    upper = np.zeros((num_points, 2))
    lower = np.zeros_like(upper)
    upper[:, 0] = x_c - thickness_distribution*np.sin(theta)
    upper[:, 1] = y_camber + thickness_distribution*np.cos(theta)
    lower[:, 0] = x_c + thickness_distribution*np.sin(theta)
    lower[:, 1] = y_camber - thickness_distribution*np.cos(theta)
    upper = np.flip(upper, axis = 0)
    coordinates = np.concatenate((upper, lower[1:]), axis = 0)

    return coordinates