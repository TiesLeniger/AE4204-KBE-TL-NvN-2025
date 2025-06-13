# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part, action
from parapy.core.validate import OneOf, Range, GE

from glidesign.geometry import lifting_surface
# Custom imports
from .lifting_surface import LiftingSurface, LiftingSection
from .fuselage import GliderFuselage
from ..analysis import ScissorPlot

class Glider(GeomBase):

    # Top-level parameters
    occupants: int = Input(1, validator = OneOf([1, 2]))                                    # Number of occupants of the glider, default value is 1
    fai_class: str = Input("std", validator = OneOf(["std", "15", "18", "20", "open"]))     # FAI class of the glider, can be "std", "15", "18", "20" or "open"
    engine_type: str = Input("No engine", validator = OneOf([                               # Glider support engine type
        "FES", "Turbo", "Self launch", "Jet", "No engine"
    ]))
    open_class_wingspan = Input(25, validator = GE(18))                                     #In case of open class glider

    # Main wing parameters
    wing_airfoil_id: str = Input()       # TODO: add validator                              # Can be NACA 4- or 5-digit or a string referencing a '.dat' file with coordinates
    wing_twist: float = Input(-2.0, validator = Range(-5.0, 5.0))                           # Twist of tip w.r.t root in [deg]
    wing_dihedral: float = Input(3.0, validator = Range(-5.0, 5.0))                         # Wing dihedral angle [deg]
    wing_sweep: float = Input(0.0, validator = Range(-10.0, 10.0))                          # Quarter chord sweep angle [deg]
    wing_incidence: float = Input(0.0, validator = Range(-5.0, 5.0))                        # Wing Incidence angle [deg]
    wing_taper: float = Input(0.5, validator = Range(0.1, 1.0))                             # Taper ratio
    flap_type: str = Input("No flaps", validator = OneOf([
        "No flaps", "Flaperon", "Discrete flap"]))
    wing_pos_long: float = Input(0.3, validator = Range(0, 1))                              #Longitudinal wing position as fraction fuselage length
    wing_pos_vert: float = Input(0.2, validator = Range(0, 1))                              #Vertical wing position as fraction max fus radius

    # Winglet parameters
    winglet_length: float = Input(0.0, validator = Range(0.0, 1.0))                         # Winglet length in [m]
    winglet_cant: float = Input(0.0, validator = Range(0.0, 90.0))                          # Cant angle of the winglet [deg] (90 deg means wing extension)
    winglet_toe: float = Input(0.0, Range(-5.0, 5.0))                                       # Toe angle of the winglet [deg]
    winglet_sweep: float = Input(0.0, Range(0.0, 30.0))                                     # Leading edge sweep of the winglet [deg]
    winglet_airfoil_id: str = Input()      # TODO: add validator                            # Airfoil profile of the winglet
    winglet_stagger: float = Input(0.0, Range(0.0, 0.8))                                    # Backwards shift of winglet LE w.r.t LE of the wingtip as fraction of the tipchord
    winglet_taper: float = Input(0.0, Range(0.1, 1.0))                                      # Winglet taper ratio

    # Horizontal tail parameters
    hor_tail_airfoil_id: float = Input()        # TODO: add validator                       # Horizontal tail airfoil profile
    hor_tail_span: float = Input(1, validator = Range(0.3, 3))                              # Horizontal tail span
    hor_tail_pos_long: float = Input(1, validator = Range(0.5, 1.3))                        #Horizontal tail position as fraction of fuselage length
    hor_tail_twist: float = Input(0, validator=Range(-5.0, 5.0))                            # Twist of tip w.r.t root in [deg]
    hor_tail_dihedral: float = Input(0.0, validator=Range(-5.0, 5.0))                       # Dihedral angle [deg]
    hor_tail_sweep: float = Input(5, validator=Range(-5.0, 5.0))                            # Quarter chord sweep angle [deg]
    hor_tail_incidence: float = Input(2, validator=Range(-5.0, 5.0))                        # Incidence angle [deg]
    hor_tail_taper: float = Input(0.5, validator=Range(0.1, 1.0))                           # Taper ratio

    # Vertical tail parameters
    ver_tail_airfoil_id: float = Input()        # TODO: add validator                       # Vertical tail airfoil profile
    ver_tail_pos_long: float = Input(1, validator=Range(0.5,1.3))                           # Vertical tail position as fraction of fuselage length
    ver_tail_height: float = Input(1.2, validator= Range(0.2, 2))                           # Height of vertical tailplane in meters
    ver_tail_sweep: float = Input(5, validator=Range(-5.0, 5.0))                            # Quarter chord sweep angle [deg]
    ver_tail_taper: float = Input(0.6, validator=Range(0.1, 1.0))                           # Taper ratio

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
            return 15 #meters
        elif self.fai_class == "15":
            return 15 #metres
        elif self.fai_class == "18":
            return 18 #metres
        elif self.fai_class == "20":
            return 20 #metres
        elif self.fai_class == "open":
            return self.open_class_wingspan                                                 # Open class has no wingspan limitation, user can define it (default = 25m)

    @Attribute
    def wing_position(self):
        return translate(self.position,
                        'x', self.wing_pos_long * self.fuselage_length,
                        'z', self.wing_pos_vert * self.fuselage_max_radius)

    @Attribute
    def wing_surface_area(self):
        return self.right_wing.area * 2

    @Attribute
    def hor_tail_surface_area(self):
        return self.right_hor_tail.area * 2

    @Attribute
    def tail_length(self):
        quarter_chord_wing = 0.25 * self.right_wing.mean_aerodynamic_chord
        quarter_chord_hor_tail = 0.25 * self.right_hor_tail.mean_aerodynamic_chord
        return (self.hor_tail_position[0] + quarter_chord_hor_tail) - (self.wing_position[0] + quarter_chord_wing)

    @Attribute
    def wing_aspect_ratio(self):
        return self.wingspan**2 / self.wing_surface_area

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
                span = self.wingspan / 2,
                twist = self.wing_twist,
                dihedral = self.wing_dihedral,
                sweep = self.wing_sweep,
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
    def right_hor_tail(self):
        return LiftingSurface(
            name = "Horizontal Tail",
            incidence_angle = 0.0,
            sections = [LiftingSection(
                idx = 0,
                root_af = "nlf1-0015",
                tip_af = "nlf1-0015",
                root_chord = 0.5,
                tip_chord = 0.3,
                span = self.hor_tail_span,
                twist = self.hor_tail_twist,
                dihedral = self.hor_tail_dihedral,
                sweep = self.hor_tail_sweep,
                sweep_loc = 0.25
            )],
            mesh_deflection = 1e-4,
            af_cst_order = 4,
            af_num_points = 40,
            af_closed_TE = True,
            position = self.hor_tail_position
        )

    @Part
    def left_hor_tail(self):
        return MirroredShape(
            shape_in = self.right_hor_tail,
            reference_point = self.hor_tail_position,
            vector1 = self.position.Vz,
            vector2 = self.position.Vx
            #mesh_deflection=self.mesh_deflection
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
                root_af = "nlf1-0015",
                tip_af = "nlf1-0015",
                root_chord = 0.8,
                tip_chord = 0.5,
                span = self.ver_tail_height,
                twist = 0.0,
                dihedral = 0.0,
                sweep = self.ver_tail_sweep,
                sweep_loc = 0.25
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

    @action
    def scissor_plot(self):
        plot = ScissorPlot(
            SM= 0.05,
            x_ac= 0.25,
            s_h= self.hor_tail_surface_area,
            s= self.wing_surface_area,
            l_h= self.tail_length,
            chord= self.right_wing.mean_aerodynamic_chord,
            cm_ac=-0.05,
            cl_a_min_h=0.7,
            cl_h=-0.2,
            cl_alpha_h=5,
            cl_alpha_a_min_h=6,
            velocity=55,
            velocity_h=55,
            wingspan = self.wingspan,
            m_tv=abs(self.hor_tail_position[-1] - self.wing_position[-1]),
            sweep_4c= np.deg2rad(self.wing_sweep),
            AR= self.wing_aspect_ratio,
        )

        # Plot
        plot.plot_scissor_plot()

        return plot

if __name__ == '__main__':
    from parapy.gui import display
    glider = Glider(name="glider")
    display(glider)