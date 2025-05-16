# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range
from unicodedata import mirrored

# Custom imports
from .lifting_surface import LiftingSurface

class Glider(GeomBase):

    # Top-level parameters
    occupants: int = Input(1, validator = OneOf([1, 2]))                                    # Number of occupants of the glider, default value is 1
    fai_class: str = Input("std", validator = OneOf(["std", "15", "18", "20", "open"]))     # FAI class of the glider, can be "std", "15", "18", "20" or "open"
    max_glide_ratio: float = Input(40.0, validator = Range(25, 75))                         # Maximum lift to drag, or glide, ratio
    min_descent_rate: float = Input(-0.3, validator = Range(-0.5, -0.2))                    # Minimum descent rate of the glider
    water_vol: float = Input(0.0, validator = Range(0.0, 250.0))                            # Water tank volume in the glider [L]
    engine_type: str = Input("No engine", validator = OneOf([                               # Glider support engine type
        "FES", "Turbo", "Self launch", "Jet", "No engine"
    ]))

    # Main wing parameters
    airfoil_id: str = Input()       # TODO: add validator                                   # Can be NACA 4- or 5-digit or a string referencing a '.dat' file with coordinates
    twist: float = Input(-2.0, validator = Range(-5.0, 0.0))                                # Twist of tip w.r.t root in [deg]
    dihedral: float = Input(0.0, validator = Range(-3.0, 3.0))                              # Wing dihedral angle [deg]
    sweep: float = Input(0.0, validator = Range(-5.0, 5.0))                                 # Quarter chord sweep angle [deg]
    taper: float = Input(0.5, validator = Range(0.1, 1.0))                                  # Taper ratio
    flap_type: str = Input("No flaps", validator = OneOf([
        "No flaps", "Flaperon", "Discrete flap"
    ]))

    # Winglet parameters
    winglet_length: float = Input(0.0, validator = Range(0.0, 1.0))                         # Winglet length in [m]
    winglet_cant: float = Input(0.0, validator = Range(0.0, 90.0))                          # Cant angle of the winglet [deg] (90 deg means wing extension)
    winglet_toe: float = Input(0.0, Range(-5.0, 5.0))                                       # Toe angle of the winglet [deg]
    winglet_sweep: float = Input(0.0, Range(0.0, 30.0))                                     # Leading edge sweep of the winglet [deg]
    winglet_airfoil_id: str = Input()      # TODO: add validator                            # Airfoil profile of the winglet
    winglet_stagger: float = Input(0.0, Range(0.0, 0.8))                                    # Backwards shift of winglet LE w.r.t LE of the wingtip as fraction of the tipchord
    winglet_taper: float = Input(0.0, Range(0.1, 1.0))                                      # Winglet taper ratio

    # Tail parameters
    hor_tail_airfoil_id: float = Input()        # TODO: add validator                       # Horizontal tail airfoil profile
    ver_tail_airfoil_id: float = Input()        # TODO: add validator                       # Vertical tail airfoil profile

    @Attribute
    def wingspan(self):
        #Define the wingspan based on FAI class limitations:
        if self.fai_class == "std":
            wingspan = 15 #meters
        elif self.fai_class == "15":
            wingspan = 15 #metres
        elif self.fai_class == "18":
            wingspan = 18 #metres
        elif self.fai_class == "20":
            wingspan = 20 #metres
        elif self.fai_class == "open":
            wingspan = Input(25) #Open class has no wingspan limitation, user can define it (default = 25m)
        return self.wingspan

    @Attribute
    def wing_position(self):
        return

    @Part
    def right_wing(self):
        return LiftingSurface(
            airfoil_id = self.airfoil_id,
            winspan = self.wingspan,
            twist = self.twist,
            dihedral = self.dihedral,
            sweep = self.sweep,
            taper=self.taper,
            flap_type = self.flap_type
        )

    @Part
    def left_wing(self):
        return MirroredShape(
            shape_in=self.right_wing,
            reference_point=self.position,
            vector1=self.position.Vz,
            vector2=self.position.Vx
            #mesh_deflection=self.mesh_deflection
        )

    @Attribute
    def hor_tail_position(self):
        return

    @Part
    def hor_tail(self):
        return LiftingSurface(
            airfoil_id = self.hor_tail_airfoil_id,
            position = self.hor_tail_position
            #TODO: add a certain default for planform parameters such that they can later be deternined in the analyses based on wing size
        )

    @Attribute
    def ver_tail_position(self):
        return

    @Part
    def vert_tail(self):
        return LiftingSurface(
            airfoil_id = self.vert_tail_airfoil_id,
            position = self.vert_tail_position
            # TODO: add a certain default for planform parameters such that they can later be deternined in the analyses based on wing size
        )

    @Attribute
    def winglet_position(self):
        return

    @Part
    def winglet(self):
        return LiftingSurface(
            position = self.winglet_position,
            span = self.winglet_length
            #TODO: add other winglet parameters
        )