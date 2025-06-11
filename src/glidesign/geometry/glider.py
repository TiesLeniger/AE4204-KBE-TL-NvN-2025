# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range
from unicodedata import mirrored

# Custom imports
from .lifting_surface import LiftingSurface, LiftingSection
from .fuselage import GliderFuselage

class Glider(GeomBase):

    # Top-level parameters
    occupants: int = Input(1, validator = OneOf([1, 2]))                                    # Number of occupants of the glider, default value is 1
    fai_class: str = Input("std", validator = OneOf(["std", "15", "18", "20", "open"]))     # FAI class of the glider, can be "std", "15", "18", "20" or "open"
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
        "No flaps", "Flaperon", "Discrete flap"]))
    wing_pos_long: float = Input(0.3, validator = Range(0, 1))                          #Longitudinal wing position as fraction fuselage length
    wing_pos_vert: float = Input(0.2, validator = Range(0, 1))                          #Vertical wing position as fraction max fus radius

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
    hor_tail_pos_long: float = Input(1, validator = Range(0.5, 1.3))             #Horizontal tail position as fraction of fuselage length

    ver_tail_airfoil_id: float = Input()        # TODO: add validator                       # Vertical tail airfoil profile
    ver_tail_pos_long: float = Input(1, validator=Range(0.5,1.3))                #Vertical tail position as fraction of fuselage length
    ver_tail_height: float = Input(1.4, validator= Range(0.2, 2))                #Height of vertical tailplane in meters

    # Fuselage parameters
    fuselage_length: float = Input(6.5, validator = Range(4.0, 12.0))
    fuselage_max_diameter: float = Input(0.7, validator = Range(0.3, 1.8))
    
    @Attribute
    def fuselage_max_radius(self):
        return self.fuselage_max_diameter/2

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
        return translate(self.position,
                        'x', self.wing_pos_long * self.fuselage_length,
                        'z', self.wing_pos_vert * self.fuselage_max_radius)

    @Part
    def right_wing(self):
        return LiftingSurface(
            name = "Main wing",
            incidence_angle = 0.0,
            sections = [LiftingSection(
                idx = 0,
                root_af = "nlf1-0015",
                tip_af = "nlf1-0015",
                root_chord = 0.5,
                tip_chord = 0.2,
                span = 7.3,
                twist = 0.0,
                dihedral = 0.0,
                sweep = 0.0,
                sweep_loc = 0.25
            )],
            mesh_deflection = 1e-4,
            af_cst_order = 5,
            af_num_points = 40,
            af_closed_TE = True,
            position = self.wing_position
        )

    @Part
    def left_wing(self):
        return MirroredShape(
            shape_in = self.right_wing,
            reference_point = self.wing_position,
            vector1 = self.position.Vz,
            vector2 = self.position.Vx
            #mesh_deflection=self.mesh_deflection
        )

    @Attribute
    def hor_tail_position(self):
        return translate(self.position,
                         'x', self.hor_tail_pos_long * self.fuselage_length,
                         'z', self.ver_tail_height)

    @Part
    def hor_tail(self):
        return LiftingSurface(
            name = "Horizontal Tail",
            incidence_angle = 0.0,
            sections = [LiftingSection(
                idx = 0,
                root_af = "naca0010",
                tip_af = "naca0010",
                root_chord = 0.4,
                tip_chord = 0.2,
                span = 1.0,
                twist = 0.0,
                dihedral = 0.0,
                sweep = 10.0,
                sweep_loc = 0.25
            )],
            mesh_deflection = 1e-4,
            af_cst_order = 5,
            af_num_points = 40,
            af_closed_TE = True,
            position = self.hor_tail_position
        )

    @Attribute
    def ver_tail_position(self):
        return rotate(translate(self.position,
                                'x', self.ver_tail_pos_long * self.fuselage_length),
                      'x', np.radians(90))

    @Part
    def vert_tail(self):
        return LiftingSurface(
            name = "Horizontal Tail",
            incidence_angle = 0.0,
            sections = [LiftingSection(
                idx = 0,
                root_af = "naca0010",
                tip_af = "naca0010",
                root_chord = 0.4,
                tip_chord = 0.2,
                span = 1.5,
                twist = 0.0,
                dihedral = 0.0,
                sweep = 10.0,
                sweep_loc = 0.0
            )],
            mesh_deflection = 1e-4,
            af_cst_order = 5,
            af_num_points = 40,
            af_closed_TE = True,
            position = self.ver_tail_position
        )

    # @Part
    # def winglet(self):
    #     return LiftingSurface(
    #         position = self.winglet_position,
    #         span = self.winglet_length
    #         #TODO: add other winglet parameters
    #     )

    @Part
    def fuselage(self):
        return GliderFuselage(
            position = self.position, #Starts at origin
            L = self.fuselage_length,
            D = self.fuselage_max_diameter,
        )

if __name__ == '__main__':
    from parapy.gui import display
    glider = Glider(name="glider")
    display(glider)