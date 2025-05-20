# Python native imports
import os

# Python third party package imports
import numpy as np

# Parapy imports
from parapy.geom import FittedCurve, Point
from parapy.core import Attribute, Part, Input
from parapy.core.validate import OneOf, Range
from kbeutils.data import airfoils

# Self-built imports
from ..core import airfoil_found, generate_naca4
from ..core.cst_curves import fit_cst_airfoil
from .ref_frame import Frame

class Airfoil(FittedCurve):

    airfoil_name: str = Input("nlf1-0015", validator = airfoil_found)
    chord: float = Input(1.0)
    cst_poly_order: int = Input(4)
    num_points: int = Input(30)
    points_spacing: str = Input('cosine', validator = OneOf(['cosine', 'constant']))
    closed_TE: bool = Input(True)

    @Attribute
    def coords(self):
        name_clean = self.airfoil_name.lower().replace(' ', '')
        if name_clean.startswith('naca') and len(name_clean) == 8:
            digits = name_clean[4:]
            try:
                m = int(digits[0])
                p = int(digits[1])
                tt = int(digits[2:])
            except ValueError:
                raise ValueError(f"Invalid format for NACA 4 series: {self.airfoil_name}")
            coordinates = generate_naca4(num_points = self.num_points, m = m, p = p, xx = tt, spacing = self.points_spacing, closed_TE = self.closed_TE)
            return [coordinates[:, 0].tolist(), coordinates[:, 1].tolist()]
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
        return fit_cst_airfoil(self.x, self.z, self.cst_poly_order)

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