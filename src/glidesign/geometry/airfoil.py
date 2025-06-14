# Python native imports
import os

# Python third party package imports
import numpy as np

# Parapy imports
from parapy.geom import FittedCurve, Point
from parapy.core import Attribute, Part, Input
from parapy.core.validate import OneOf, Range, GreaterThan
from parapy.gui import display
from kbeutils.data import airfoils
from kbeutils.geom import Naca4AirfoilCurve, Naca5AirfoilCurve

# Self-built imports
from ..core import airfoil_found, NACA4_PATTERN, NACA5_PATTERN
from ..core.cst_curves import fit_cst_airfoil
from .ref_frame import Frame

class Airfoil(FittedCurve):

    airfoil_name: str = Input("nlf1-0015", validator = airfoil_found)
    chord: float = Input(1.0, validator = GreaterThan(0.0))
    cst_poly_order: int = Input(5, validator = GreaterThan(0))
    num_points: int = Input(200, validator = GreaterThan(0))
    mesh_deflection: float = Input(1e-5, validator = GreaterThan(0.0))

    @Attribute
    def coords(self):
        if NACA4_PATTERN.match(self.airfoil_name):
            digits = NACA4_PATTERN.match(self.airfoil_name).group(1)
            naca4curve = Naca4AirfoilCurve(designation = digits, mesh_deflection = 1e-6, hidden = True)
            coordinates = np.array(np.array(naca4curve.coordinates))
            assert coordinates.shape == (naca4curve.n_points, 3), f"Expected shape ({naca4curve.n_points}, 3), got {coordinates.shape}"
            coordinates = np.round(coordinates / naca4curve.chord_length, decimals = 5)
            return [coordinates[:, 0].tolist(), coordinates[:, 2].tolist()]
        elif NACA5_PATTERN.match(self.airfoil_name):
            digits = NACA5_PATTERN.match(self.airfoil_name).group(1)
            naca5curve = Naca5AirfoilCurve(designation = digits, mesh_deflection = 1e-6, hidden = True)
            coordinates = np.array(np.array(naca5curve.coordinates))
            assert coordinates.shape == (naca5curve.n_points, 3), f"Expected shape ({naca5curve.n_points}, 3), got {coordinates.shape}"
            coordinates = np.round(coordinates / naca5curve.chord_length, decimals = 5)
            return [coordinates[:, 0].tolist(), coordinates[:, 2].tolist()]
        else:
            in_input = os.path.exists(os.path.join(os.getcwd(), "input", "airfoils", self.airfoil_name + ".dat"))
            path_to_af_file = os.path.join(os.getcwd(), "input", "airfoils", self.airfoil_name + ".dat") if in_input else os.path.join(
                airfoils.__path__[0], self.airfoil_name + ".dat")
            with open(path_to_af_file, 'r') as file:
                x_coords_lst = []
                z_coords_lst = []
                skipped_lines = 0
                for line in file:
                    try:
                        x, z = line.split(maxsplit = 1)
                        x_coords_lst.append(float(x.strip()))
                        z_coords_lst.append(float(z.strip()))
                    except ValueError:
                        if skipped_lines == 0:
                            skipped_lines = 1
                            continue
                        else:
                            break
            return [x_coords_lst, z_coords_lst]

    @Attribute
    def points(self) -> list[Point]:
        return [self.position.translate("x", x * self.chord,
                                               "z", z * self.chord).location
                        for x, z in zip(self.coords[0], self.coords[1])]
    
    @Attribute
    def cst_coefficients(self):
        return fit_cst_airfoil(np.array(self.coords[0]), np.array(self.coords[1]), self.cst_poly_order)

    @Attribute
    def cst_coeff_u(self):
        return self.cst_coefficients[0]
    
    @Attribute
    def cst_coeff_l(self):
        return self.cst_coefficients[1]
    
    @Attribute
    def z_te_u(self):
        return self.cst_coefficients[2]
    
    @Attribute
    def z_te_l(self):
        return self.cst_coefficients[3]
    
    @Part
    def frame(self):
        return Frame(pos = self.position, hidden = False)