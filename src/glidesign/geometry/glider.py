# Python native imports
import warnings
from wsgiref.validate import validator

# Python third party imports
import numpy as np
import matlab

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part, action
from parapy.core.widgets import Dropdown
from parapy.core.validate import OneOf, Range, GE, Validator, GreaterThan
import kbeutils.avl as avl

# Custom imports
from .lifting_surface import LiftingSurface, LiftingSection
from .fuselage import GliderFuselage

from ..analysis import ScissorPlot
from ..core import airfoil_found, convert_matlab_dict
from ..external import MATLAB_Q3D_ENGINE, Q3DData
from ..analysis import ScissorPlot, WeightAndBalance
from ..core import airfoil_found


class Glider(GeomBase):

    # Top-level parameters
    occupants: int = Input(1, widget = Dropdown([1, 2]), validator = OneOf([1, 2]))                         # Number of occupants of the glider, default value is 1
    fai_class: str = Input(
        "standard class",                                                                                   # FAI class of the glider, can be "std", "15", "18", "20" or "open" 
        widget = Dropdown(["standard class", "15m class", "18m class", "20m class", "open class"]),
        validator = OneOf(["standard class", "15m class", "18m class", "20m class", "open class"]))     
    open_class_wingspan = Input(25, validator = GE(18))                                                     #In case of open class glider
    cl_cr: float = Input(0.3, validator = GreaterThan(0.0))

    min_pilot_mass = Input(50, validator= Range(50, 90))                 #Minimum allowed mass for pilot
    max_pilot_mass = Input(100, validator= Range(80, 120))               #Maximum allowed pilot mass
    glider_structure_material = Input("Carbon fibre", widget = Dropdown(["Carbon fibre", "Glass fibre"]))

    # Main wing parameters
    wing_airfoil_id: str = Input('nlf1-0015', validator = airfoil_found)                    # Can be NACA 4- or 5-digit or a string referencing a '.dat' file with coordinates
    wing_twist: float = Input(-2.0, validator = Range(-5.0, 5.0))                           # Twist of tip w.r.t root in [deg]
    wing_dihedral: float = Input(3.0, validator = Range(-5.0, 5.0))                         # Wing dihedral angle [deg]
    wing_sweep: float = Input(0.0, validator = Range(-10.0, 10.0))                          # Quarter chord sweep angle [deg]
    wing_incidence: float = Input(0.0, validator = Range(-5.0, 5.0))                        # Wing Incidence angle [deg]
    wing_taper: float = Input(0.5, validator = Range(0.1, 1.0))                             # Taper ratio
    wing_pos_long: float = Input(0.25, validator = Range(0.0, 1.0))
    wing_pos_vert: float = Input(0.2, validator = Range(0.0, 1.0))

    # Winglet parameters
    has_winglet: bool = Input(True)                                                         # Boolean for adding a winglet
    winglet_length: float = Input(0.3, validator = Range(0.0, 1.0))                         # Winglet length in [m]
    winglet_cant: float = Input(5.0, validator = Range(0.0, 90.0))                          # Cant angle of the winglet [deg] (90 deg means wing extension)
    winglet_toe: float = Input(0.5, Range(-5.0, 5.0))                                       # Toe angle of the winglet [deg]
    winglet_sweep: float = Input(30.0, Range(0.0, 30.0))                                    # Leading edge sweep of the winglet [deg]
    winglet_tip_af: str = Input('NACA 0010', validator = airfoil_found)                     # Airfoil profile of the winglet
    winglet_taper: float = Input(0.5, Range(0.1, 1.0))                                      # Winglet taper ratio

    # Horizontal tail parameters
    hor_tail_airfoil_id: float = Input('NACA 0010', validator = airfoil_found)              # Horizontal tail airfoil profile
    hor_tail_span: float = Input(2.5, validator = Range(0.3, 15))                              # Horizontal tail span

    hor_tail_pos_long: float = Input(1, validator = Range(0.5, 1.3))                        # Horizontal tail position as fraction of fuselage length
    hor_tail_twist: float = Input(0, validator=Range(-5.0, 5.0))                            # Twist of tip w.r.t root in [deg]
    hor_tail_dihedral: float = Input(0.0, validator=Range(-5.0, 5.0))                       # Dihedral angle [deg]
    hor_tail_sweep: float = Input(5, validator=Range(-5.0, 5.0))                            # Quarter chord sweep angle [deg]
    hor_tail_incidence: float = Input(2, validator=Range(-5.0, 5.0))                        # Incidence angle [deg]
    hor_tail_taper: float = Input(0.55, validator=Range(0.1, 1.0))                           # Taper ratio
    SM: float = Input(0.1, validator = Range(0, 0.2))                                       # Stability margin

    # Vertical tail parameters
    ver_tail_airfoil_id: float = Input('NACA 0012', validator = airfoil_found)              # Vertical tail airfoil profile
    ver_tail_root_chord: float = Input(0.9, validator = GreaterThan(0.0))                   # Vertical tail root chord
    ver_tail_pos_long: float = Input(1, validator=Range(0.5,1.3))                           # Vertical tail position as fraction of fuselage length
    ver_tail_height: float = Input(1.2, validator= Range(0.2, 2))                           # Height of vertical tailplane in meters
    ver_tail_sweep: float = Input(5, validator=Range(-5.0, 5.0))                            # Quarter chord sweep angle [deg]
    ver_tail_taper: float = Input(0.6, validator=Range(0.1, 1.0))                           # Taper ratio

    @Input
    def current_pilot_mass(self):
        #Return the average for current pilot mass estimate
        #Also adjustable
        return (self.max_pilot_mass + self.min_pilot_mass) / 2

    # Fuselage parameters
    @Input(validator = Range(0.3, 1.8))
    def fuselage_max_diameter(self):
        val = 0.0
        if self.occupants == 1:
            val = 0.7
        if self.occupants == 2:
            val = 1.0
        return val

    @Input(validator = Range(4.0, 10.0))
    def fuselage_length(self):
        val = 0.0
        if self.occupants == 1:
            val = 6.5
        if self.occupants == 2:
            val = 8.5
        return val

    #Glider Attributes and Parts

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
    def hor_tail_span(self):
        #Define the tailspan based on wingspan:
        if self.wingspan <= 18:
            return 2.5 #meters
        else:
            return 3.0 #metres

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
        return self.right_wing.half_area * 2

    @Attribute
    def hor_tail_surface_area(self):
        return self.right_hor_tail.half_area * 2

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
            root_af = self.wing_airfoil_id,
            tip_af = self.wing_airfoil_id,
            root_chord = 0.9,
            taper = 0.4,
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
            has_winglet = self.has_winglet,
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
    
    @Attribute(in_tree = True)
    def main_wing_avl_surface(self):
        return avl.Surface(
            name = "Main Wing",
            n_chordwise = Input(12, validator = GreaterThan(0)),
            chord_spacing = avl.Spacing.cosine,
            n_spanwise = Input(24, validator = GreaterThan(0)),
            span_spacing = avl.Spacing.cosine,
            y_duplicate = self.right_wing.position.point[1],
            sections = [profile.avl_section for profile in self.right_wing.profiles]   
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
            taper = 0.55,
            span = self.hor_tail_span / 2,
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
            has_winglet = False,
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
    
    @Attribute(in_tree = True)
    def horizontal_tail_avl_surface(self):
        return avl.Surface(
            name = "Horizontal tail",
            n_chordwise = Input(8, validator = GreaterThan(0)),
            chord_spacing = avl.Spacing.cosine,
            n_spanwise = Input(16, validator = GreaterThan(0)),
            span_spacing = avl.Spacing.cosine,
            y_duplicate = self.right_hor_tail.position.point[1],
            sections = [profile.avl_section for profile in self.right_hor_tail.profiles]   
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
            root_af = "NACA 0014",
            tip_af = "NACA 0014",
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
            has_winglet = False,
            winglet_tip_af = 'NACA 0010'
        )
    
    @Attribute(in_tree = True)
    def vertical_tail_avl_surface(self):
        return avl.Surface(
            name = "Vertical tail",
            n_chordwise = Input(8, validator = GreaterThan(0)),
            chord_spacing = avl.Spacing.cosine,
            n_spanwise = Input(12, validator = GreaterThan(0)),
            span_spacing = avl.Spacing.cosine,
            y_duplicate = None,
            sections = [profile.avl_section for profile in self.vert_tail.profiles]   
        )

    @Part
    def fuselage(self):
        return GliderFuselage(
            position = self.position, #Starts at origin
            L = self.fuselage_length,
            D = self.fuselage_max_diameter,
            mesh_deflection=1e-4,
            wing_position = self.wing_position,
        )

    @Attribute
    def glider_x_cog(self):
        glider_wb = WeightAndBalance(min_pilot_mass=self.min_pilot_mass,
                                     max_pilot_mass=self.max_pilot_mass,
                                     current_pilot_mass=self.current_pilot_mass,
                                     glider_structure_material=self.glider_structure_material,
                                     right_wing_cog=self.right_wing.cog,
                                     left_wing_cog=self.right_wing.cog,
                                     vertical_tail_cog= self.vert_tail.cog,
                                     right_hor_tail_cog=self.right_hor_tail.cog,
                                     left_hor_tail_cog=self.left_hor_tail.cog,
                                     fuselage_cog=self.fuselage.fuselage_solid.cog,
                                     pilot_cog=self.fuselage.canopy_ellipse.center,
                                     right_wing_volume=self.right_wing.volume,
                                     left_wing_volume = self.right_wing.volume, #Mirrored has no volume for some reason (assume symmetrical)
                                     vertical_tail_volume = self.vert_tail.volume,
                                     right_hor_tail_volume = self.right_hor_tail.volume,
                                     left_hor_tail_volume = self.right_hor_tail.volume, #Mirrored has no volume for some reason (assume symmetrical)
                                     fuselage_volume = self.fuselage.fuselage_solid.volume,
        )
        current_x_cog = glider_wb.get_current_cog
        fwd_limit_x_cog = glider_wb.get_fwd_limit_cog
        bwd_limit_x_cog = glider_wb.get_bwd_limit_cog
        return current_x_cog, fwd_limit_x_cog, bwd_limit_x_cog

    @Attribute
    def glider_empty_mass(self):
        glider_wb = WeightAndBalance(min_pilot_mass=self.min_pilot_mass,
                                     max_pilot_mass=self.max_pilot_mass,
                                     current_pilot_mass=self.current_pilot_mass,
                                     glider_structure_material=self.glider_structure_material,
                                     right_wing_cog=self.right_wing.cog,
                                     left_wing_cog=self.right_wing.cog,
                                     vertical_tail_cog= self.vert_tail.cog,
                                     right_hor_tail_cog=self.right_hor_tail.cog,
                                     left_hor_tail_cog=self.left_hor_tail.cog,
                                     fuselage_cog=self.fuselage.fuselage_solid.cog,
                                     pilot_cog=self.fuselage.canopy_ellipse.center,
                                     right_wing_volume=self.right_wing.volume,
                                     left_wing_volume = self.right_wing.volume, #Mirrored has no volume for some reason (assume symmetrical)
                                     vertical_tail_volume = self.vert_tail.volume,
                                     right_hor_tail_volume = self.right_hor_tail.volume,
                                     left_hor_tail_volume = self.right_hor_tail.volume, #Mirrored has no volume for some reason (assume symmetrical)
                                     fuselage_volume = self.fuselage.fuselage_solid.volume
        )
        current_empty_mass = glider_wb.get_total_empty_mass
        return current_empty_mass

    @Attribute
    def glider_list_of_masses(self):
        glider_wb = WeightAndBalance(min_pilot_mass=self.min_pilot_mass,
                                     max_pilot_mass=self.max_pilot_mass,
                                     current_pilot_mass=self.current_pilot_mass,
                                     glider_structure_material=self.glider_structure_material,
                                     right_wing_cog=self.right_wing.cog,
                                     left_wing_cog=self.right_wing.cog,
                                     vertical_tail_cog= self.vert_tail.cog,
                                     right_hor_tail_cog=self.right_hor_tail.cog,
                                     left_hor_tail_cog=self.left_hor_tail.cog,
                                     fuselage_cog=self.fuselage.fuselage_solid.cog,
                                     pilot_cog=self.fuselage.canopy_ellipse.center,
                                     right_wing_volume=self.right_wing.volume,
                                     left_wing_volume = self.right_wing.volume, #Mirrored has no volume for some reason (assume symmetrical)
                                     vertical_tail_volume = self.vert_tail.volume,
                                     right_hor_tail_volume = self.right_hor_tail.volume,
                                     left_hor_tail_volume = self.right_hor_tail.volume, #Mirrored has no volume for some reason (assume symmetrical)
                                     fuselage_volume = self.fuselage.fuselage_solid.volume
        )
        list_of_masses = glider_wb.list_of_masses
        return list_of_masses

    @action(button_label = "plot")
    def scissor_plot(self):
        plot = ScissorPlot(
            wing_x_location = self.wing_position[0], #Should technically be MAC position
            current_x_cog= self.glider_x_cog[0],
            fwd_limit_x_cog= self.glider_x_cog[1],
            bwd_limit_x_cog= self.glider_x_cog[2],
            SM= self.SM,
            x_ac= 0.25,
            s= self.wing_surface_area,
            l_h= self.tail_length,
            chord= self.right_wing.mean_aerodynamic_chord,
            cm_ac=-0.05,
            cl_a_min_h=0.7,
            cl_h=-0.2,
            cl_alpha_h=5,
            cl_alpha_a_min_h=6,
            velocity=self.Q3D_params.velocity,
            velocity_h=self.Q3D_params.velocity,
            wingspan = self.wingspan,
            m_tv=abs(self.hor_tail_position[-1] - self.wing_position[-1]),
            sweep_4c= np.deg2rad(self.wing_sweep),
            AR= self.wing_aspect_ratio,
            hor_tail_span= self.hor_tail_span,
            hor_tail_taper = self.hor_tail_taper
        )
        self.right_hor_tail.root_chord = plot.c_h_root_auto_scaled

        # Plot
        plot.plot_scissor_plot()

        return plot
    
    @Part
    def Q3D_params(self):
        return Q3DData()

    @action(label = "Run Q3D")
    def q3d_data(self):
        """All inputs and results from running Q3D (MATLAB)"""
        self.q3d_res = MATLAB_Q3D_ENGINE.run_q3d_cst(
            self.right_wing.q3d_planform_geom,
            self.right_wing.q3d_cst_airfoils,
            self.right_wing.q3d_eta_airfoils,
            matlab.double(self.right_wing.incidence_angle),
            self.Q3D_params.mach_number,
            self.Q3D_params.reynolds_number,
            self.Q3D_params.velocity,
            self.Q3D_params.alpha,
            self.Q3D_params.altitude,
            self.Q3D_params.density
        )

    @Attribute
    def q3d_wing_data(self):
        result = getattr(self, "q3d_res", None)
        if result is not None:
            return convert_matlab_dict(result["Wing"])
        else:
            return "Evaluate aerodynamics to view property"

    @Attribute
    def q3d_section_data(self):
        result = getattr(self, "q3d_res", None)
        if result is not None:
            return convert_matlab_dict(result["Section"])
        else:
            return "Evaluate aerodynamics to view property"

    @Attribute
    def wing_cl(self) -> float:
        result = getattr(self, "q3d_res", None)
        if result is not None:
            return result["CLwing"]
        else:
            return "Evaluate aerodynamics to view property"

    @Attribute
    def wing_cd(self) -> float:
        result = getattr(self, "q3d_res", None)
        if result is not None:
            return result["CDwing"]
        else:
            return "Evaluate aerodynamics to view property"

    @Attribute
    def wing_cm(self) -> float:
        result = getattr(self, "q3d_res", None)
        if result is not None:
            return result["CMwing"]
        else:
            return "Evaluate aerodynamics to view property"
        
    @Attribute(in_tree = True)
    def avl_surfaces(self):
        return [self.main_wing_avl_surface, self.horizontal_tail_avl_surface, self.vertical_tail_avl_surface]
    
    @Part
    def avl_configuration(self):
        """Configurations are made separately for each Mach number that is provided."""
        return avl.Configuration(name='cruise analysis',
                                 reference_area=self.wing_surface_area,
                                 reference_span=self.wingspan,
                                 reference_chord=self.right_wing.mean_aerodynamic_chord,
                                 reference_point=self.position.point,
                                 surfaces=self.avl_surfaces,
                                 mach= 0.0)
    
    @Attribute
    def avl_settings(self):
        return {'alpha': avl.Parameter(name='alpha',
                                         setting='CL',
                                         value=self.cl_cr)}
    
    @Part
    def avl_case(self):
        """avl case definition using the avl_settings dictionary defined above"""
        return avl.Case(name='fixed_cl',  # name _must_ correspond to type of case
                        settings=self.avl_settings)
    
    @Part
    def avl_analysis(self):
        return avl.Interface(configuration=self.avl_configuration,
                             # note: AVL always expects a list of cases!
                             cases=[self.avl_case])
    
    @Attribute
    def full_aircraft_moment_coefficient(self):
        return {result['Name']: result['Totals']['CMtot']
                for case_name, result in self.avl_analysis.results.items()}

if __name__ == '__main__':
    from parapy.gui import display
    glider = Glider(name="glider")
    display(glider)