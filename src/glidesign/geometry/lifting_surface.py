# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range, GreaterThan

# Custom imports
from ..core import airfoil_found
from .airfoil import Airfoil

class LiftingSection(GeomBase):

    root_airfoil_id: str = Input(validator = airfoil_found)             # Inner airfoil name
    tip_airfoil_id: str = Input(validator = airfoil_found)              # Outer airfoil name
    root_chord: float = Input(validator = GreaterThan(0))               # Root chord [m]
    taper_ratio: float = Input(validator = Range(0.0, 1.0))             # Taper ratio of the section
    span: float = Input(15.0, validator = GreaterThan(0))               # Section span [m]
    twist: float = Input(0.0, validator = Range(-5.0, 5.0))             # Section twist (tip w.r.t root) [deg]
    dihedral: float = Input(0.0, validator = Range(-3.0, 5.0))          # Section dihedral [deg]
    sweep: float = Input(0.0, validator = Range(-10.0, 20.0))           # Section sweep [deg]
    sweep_loc: float = Input(0.25, validator = Range(0.0, 1.0))         # Chord normalized sweep location of the section (0 is LE, 1 is TE)
    control_surface: str = Input("None", OneOf(                         # Control surface on the section
        "None", "Aileron", "Flaperon", "Flap", "Airbrake"))
    
    @Attribute
    def tip_chord(self):
        return self.root_chord * self.taper_ratio

    @Attribute
    def root_position(self):
        pass

    @Part
    def airfoil_in(self):
        return Airfoil(airfoil_name = self.in_airfoil_id,
                       chord = self.root_chord,
                       )

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