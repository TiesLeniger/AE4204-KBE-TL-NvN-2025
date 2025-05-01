# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range

class Glider(GeomBase):

    # Top-level parameters
    occupants: int = Input(1, validator = OneOf([1, 2]))                                    # Number of occupants of the glider, default value is 1
    fai_class: str = Input("std", validator = OneOf(["std", "15", "18", "20", "open"]))     # FAI class of the glider, can be "std", "15", "18", "20" or "open"
    max_glide_ratio: float = Input(40.0, validator = Range(25, 75))                         # Maximum lift to drag, or glide, ratio
    min_descent_rate: float = Input(-0.3, validator = Range(-0.5, -0.2))                    # Minimum descent rate of the glider
    water_vol: float = Input(0.0, validator = Range(0.0, 250.0))                            # Water tank volume in the glider [L]

    # Main wing parameters
    airfoil_id: str = Input()       # TODO: add validator                                   # Can be NACA 4- or 5-digit or a string referencing a '.dat' file with coordinates
    twist: float = Input(-2.0, validator = Range(-5, 0))                                    # Twist of tip w.r.t root in [deg]
    