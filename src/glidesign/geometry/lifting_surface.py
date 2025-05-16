# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range

class LiftingSurface(GeomBase):
    
    name: str = Input()

    airfoil_id: float = input()
    span: float = input()
    twist: float = input()
    dihedral: float = input()
    sweep: float  = input()
    taper: float = input()
    flap_type: str = input()

    root_airfoil_id: float = Input()
    tip_airfoil_id: float = Input()

    @Attribute
    def root_position(self):
        #....
        return np.array([0, 0, 0])

    @Attribute
    def tip_position(self):
        #.....
        return np.array([0, 0, 0])

    @Attribute
    def wing_area(self):
        return

    @Attribute
    def aspect_ratio(self):
        return

