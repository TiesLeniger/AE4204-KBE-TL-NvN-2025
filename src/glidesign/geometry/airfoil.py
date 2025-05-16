# Python native imports
import os

# Python third party package imports

# Parapy imports
from parapy.geom import FittedCurve, Point
from parapy.core import Attribute, Part, Input
from parapy.core.validate import OneOf, Range
from kbeutils.data import airfoils

# Self-built imports
from ..core import airfoil_found

class Airfoil(FittedCurve):

    airfoil_name: str = Input("nlf1-0015", validator = airfoil_found)
    chord: float = Input(1.0)
    thickness_ratio: float = Input(0.15, validator = Range(0.08, 0.2))
    cst_poly_order: int = Input(4)

    @Attribute
    def coords(self):

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
                        skipped_lines == 1
                        continue
                    else:
                        break
        return [x_coords_lst, z_coords_lst]

    @Attribute
    def points(self) -> list[Point]:
        points_list = [self.position.translate("x", x * self.chord,
                                               "z", z * self.chord).location
                        for x, z in zip(self.coords[0], self.coords[1])]
        return points_list 
    

