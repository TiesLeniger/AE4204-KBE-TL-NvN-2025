# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part, action
from parapy.core.widgets import Dropdown
from parapy.core.validate import OneOf, Range, GE, GreaterThan

from glidesign.geometry import lifting_surface
# Custom imports
from .lifting_surface import LiftingSurface, LiftingSection
from .fuselage import GliderFuselage
from ..analysis import ScissorPlot
from ..core import airfoil_found

class Glider(GeomBase):

    # Top-level parameters
    occupants: int = Input(1, validator = OneOf([1, 2]))                                    # Number of occupants of the glider, default value is 1
    fai_class: str = Input("standard class", widget = Dropdown(["standard class", "15m class", "18m class", "20m class", "open class"]) ,validator = OneOf(["standard class", "15m class", "18m class", "20m class", "open class"]))     # FAI class of the glider, can be "std", "15", "18", "20" or "open"
    open_class_wingspan = Input(25, validator = GE(18))                                     #In case of open class glider

    # Main wing parameters
    wing_airfoil_id: str = Input('NACA 0010', validator = airfoil_found)                    # Can be NACA 4- or 5-digit or a string referencing a '.dat' file with coordinates
    wing_twist: float = Input(-2.0, validator = Range(-5.0, 5.0))                           # Twist of tip w.r.t root in [deg]
    wing_dihedral: float = Input(3.0, validator = Range(-5.0, 5.0))                         # Wing dihedral angle [deg]
    wing_sweep: float = Input(0.0, validator = Range(-10.0, 10.0))                          # Quarter chord sweep angle [deg]
    wing_incidence: float = Input(0.0, validator = Range(-5.0, 5.0))                        # Wing Incidence angle [deg]
    wing_taper: float = Input(0.5, validator = Range(0.1, 1.0))                             # Taper ratio
    wing_pos_long: float = Input(0.3, validator = Range(0.0, 1.0))
    wing_pos_vert: float = Input(0.2, validator = Range(0.0, 1.0))

    # Winglet parameters
    winglet: bool = Input(True)                                                            # Boolean for adding a winglet
    winglet_length: float = Input(0.3, validator = Range(0.0, 1.0))                         # Winglet length in [m]
    winglet_cant: float = Input(5.0, validator = Range(0.0, 90.0))                          # Cant angle of the winglet [deg] (90 deg means wing extension)
    winglet_toe: float = Input(0.5, Range(-5.0, 5.0))                                       # Toe angle of the winglet [deg]
    winglet_sweep: float = Input(30.0, Range(0.0, 30.0))                                     # Leading edge sweep of the winglet [deg]
    winglet_tip_af: str = Input('NACA 0010', validator = airfoil_found)                     # Airfoil profile of the winglet
    winglet_taper: float = Input(0.5, Range(0.1, 1.0))                                      # Winglet taper ratio

    # Horizontal tail parameters
    hor_tail_airfoil_id: float = Input('NACA 0010', validator = airfoil_found)              # Horizontal tail airfoil profile
    hor_tail_span: float = Input(1, validator = Range(0.3, 3))                              # Horizontal tail span
    hor_tail_pos_long: float = Input(1, validator = Range(0.5, 1.3))                        # Horizontal tail position as fraction of fuselage length
    hor_tail_twist: float = Input(0, validator=Range(-5.0, 5.0))                            # Twist of tip w.r.t root in [deg]
    hor_tail_dihedral: float = Input(0.0, validator=Range(-5.0, 5.0))                       # Dihedral angle [deg]
    hor_tail_sweep: float = Input(5, validator=Range(-5.0, 5.0))                            # Quarter chord sweep angle [deg]
    hor_tail_incidence: float = Input(2, validator=Range(-5.0, 5.0))                        # Incidence angle [deg]
    hor_tail_taper: float = Input(0.5, validator=Range(0.1, 1.0))                           # Taper ratio
    SM: float = Input(0.1, validator = Range(0, 0.2))                                       # Stability margin

    # Vertical tail parameters
    ver_tail_airfoil_id: float = Input(validator = airfoil_found)                           # Vertical tail airfoil profile
    ver_tail_root_chord: float = Input(0.5, validator = GreaterThan(0.0))                   # Vertical tail root chord
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
        if self.fai_class == "standard class":
            return 15 #meters
        elif self.fai_class == "15m class":
            return 15 #metres
        elif self.fai_class == "18m class":
            return 18 #metres
        elif self.fai_class == "20m class":
            return 20 #metres
        elif self.fai_class == "open class":
            return self.open_class_wingspan                                                 # Open class has no wingspan limitation, user can define it (default = 25m)

    @Attribute
    def max_to_mass(self):
        #Define the MTOM based on FAI class limitations:
        if self.fai_class == "standard class":
            return 525 #kg
        elif self.fai_class == "15m class":
            return 525 #kg
        elif self.fai_class == "18m class":
            return 600 #kg
        elif self.fai_class == "20m class":
            return 800 #kg
        elif self.fai_class == "open class":
            return 850 #kg

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
            root_af = "nlf1-0015",
            tip_af = "nlf1-0015",
            root_chord = 0.5,
            taper = 0.5,
            span = 7.5,
            twist = -2.0,
            sweep = 2.0,
            sweep_loc = 0.25,
            dihedral = 2.0,
            incidence_angle = 1.0,
            num_sections = 1,
            mesh_deflection = 1e-4,
            af_cst_order = 5,
            af_num_points = 40,
            af_closed_TE = True,
            position = self.wing_position,
            winglet = self.winglet,
            winglet_length = self.winglet_length,
            winglet_cant = self.winglet_cant,
            winglet_toe = self.winglet_toe,
            winglet_sweep = self.winglet_sweep,
            winglet_tip_af = self.winglet_tip_af,
            winglet_taper = self.winglet_taper
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
            name = "Horizontal tail",
            root_af = self.hor_tail_airfoil_id,
            tip_af = self.hor_tail_airfoil_id,
            root_chord = self.vert_tail.profiles[-1].chord + 0.1,
            taper = self.hor_tail_taper,
            span = self.hor_tail_span,
            twist = self.hor_tail_twist,
            sweep = self.hor_tail_sweep,
            sweep_loc = 0.25,
            dihedral = self.hor_tail_dihedral,
            incidence_angle = self.hor_tail_incidence,
            num_sections = 1,
            mesh_deflection = 1e-4,
            af_cst_order = 5,
            af_num_points = 40,
            af_closed_TE = True,
            position = self.hor_tail_position,
            winglet = False,
            winglet_tip_af = 'NACA 0010'
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
            name = "Vertical tail",
            root_af = "NACA 0012",
            tip_af = "NACA 0012",
            root_chord = self.ver_tail_root_chord,
            taper = self.ver_tail_taper,
            span = self.ver_tail_height,
            twist = 0.0,
            sweep = self.ver_tail_sweep,
            sweep_loc = 0.0,
            dihedral = 0.0,
            incidence_angle = 0.0,
            num_sections = 1,
            mesh_deflection = 1e-4,
            af_cst_order = 5,
            af_num_points = 40,
            af_closed_TE = True,
            position = self.ver_tail_position,
            winglet = False,
            winglet_tip_af = 'NACA 0010'
        )

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
            SM= self.SM,
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