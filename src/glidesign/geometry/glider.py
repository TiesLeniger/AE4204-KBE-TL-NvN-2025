# Python native imports
from pathlib import Path

# Python third party imports
import numpy as np
import matlab

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape, Point
from parapy.core import Input, Attribute, Part, action, child
from parapy.core.widgets import Dropdown
from parapy.core.validate import OneOf, Range, GE, Validator, GreaterThan
from parapy.exchange.step import STEPWriter
import kbeutils.avl as avl

# Custom imports
from .lifting_surface import LiftingSurface, LiftingSection
from .fuselage import GliderFuselage

from ..analysis import ScissorPlot
from ..core import airfoil_found, convert_matlab_dict
from ..external import MATLAB_Q3D_ENGINE, Q3DData
from ..analysis import ScissorPlot, WeightAndBalance, size_tail
from ..core import airfoil_found, G0, RHO0

# Constants
CS22_EMPTY_STALL_SPEED = 90/3.6                         # CS22 specified stall speed with airbrakes retracted, water balast empty (pilot + empty weight), [m/s]
CL_MAX_ESTIMATE = 1.3                                   # [-], estimate of CL max that a high aspect glider wing can generate

class Glider(GeomBase):
    # Top-level parameters
    occupants: int = Input(1, widget = Dropdown([1, 2]), validator = OneOf([1, 2]))                         # Number of occupants of the glider, default value is 1
    fai_class: str = Input(
        "standard class",                                                                                   # FAI class of the glider, can be "std", "15", "18", "20" or "open" 
        widget = Dropdown(["standard class", "15m class", "18m class", "20m class", "open class"]),
        validator = OneOf(["standard class", "15m class", "18m class", "20m class", "open class"]))     

    open_class_wingspan = Input(25, validator = GE(18))                                                     # In case of open class glider
    min_pilot_mass = Input(70, validator= Range(60, 90))                                                    # Minimum allowed mass for pilot
    max_pilot_mass = Input(110, validator= Range(80, 120))                                                  # Maximum allowed pilot mass
    glider_structure_material = Input("Carbon fibre", widget = Dropdown(["Carbon fibre", "Glass fibre"]))

    wing_taper: float = Input(0.4, validator = Range(0.0, 1.0, incl_min = False))
    wing_pos_long: float = Input(0.28, validator = Range(0.0, 1.0))
    wing_pos_vert: float = Input(0.22, validator = Range(0.0, 1.0))
    wing_avl_n_chordwise = Input(12, validator = GreaterThan(0))                         # Vertical tail chordwise elements for avl
    wing_avl_n_spanwise = Input(24, validator = GreaterThan(0))                         # Vertical tail spanwise elements for avl

    @Input(validator = GreaterThan(80.0/3.6))
    def cruise_speed(self):
        #Define the MTOM based on FAI class limitations:
        if self.fai_class == "standard class":
            cruise_speed = 140                          # [km/h]
        elif self.fai_class == "15m class": 
            cruise_speed = 160
        elif self.fai_class == "18m class":
            cruise_speed = 180
        elif self.fai_class == "20m class":
            cruise_speed = 160
        elif self.fai_class == "open class":
            cruise_speed = 180
        return cruise_speed/3.6

    @Input(validator = Range(60.0, 120.0))
    def current_pilot_mass(self) -> float:
        return (self.max_pilot_mass + self.min_pilot_mass) / 2                              # Return average of min and max pilot mass

    # Inputs used for main wing sizing
    @Input
    def wing_loading(self) -> float:
        return (1/2) * RHO0 * CS22_EMPTY_STALL_SPEED**2 * CL_MAX_ESTIMATE * (1/G0)
    
    @Input(validator = GE(0.0))
    def wing_root_chord(self) -> float:
        return (2*self.wing_surface_area)/(self.wing_span*(1 + self.wing_taper))

    @Input
    def has_winglet(self):
        return False if self.fai_class == "standard class" else True

    # Horizontal tail parameters
    hor_tail_overhang: float = Input(0.0, validator = GE(0.0))                              # distance in x between LE of root of horizontal tail and LE of tip of vertical tail
    Sh_S: float = Input(0.11, validator = Range(0.0, 1.0, incl_min=False))                  # Ratio of horizotal tail surface area to main wing surface area
    hor_tail_taper: float = Input(0.6, validator=Range(0.1, 1.0))                           # Taper ratio
    hor_tail_avl_n_chordwise = Input(8, validator = GreaterThan(0))                         # Vertical tail chordwise elements for avl
    hor_tail_avl_n_spanwise = Input(16, validator = GreaterThan(0))                         # Vertical tail spanwise elements for avl

    SM: float = Input(0.1, validator = Range(0, 0.2))                                       # Stability margin
    
    @Input
    def hor_tail_span(self):
        if self.wing_span <= 18:
            return 2.5
        else:
            return 3.0
    
    @Input
    def hor_tail_surface_area(self):
        return self.Sh_S * self.wing_surface_area

    @Input
    def hor_tail_root_chord(self):
        return (2*self.hor_tail_surface_area)/(self.hor_tail_span * (1 + self.hor_tail_taper))

    # Vertical tail parameters
    ver_tail_aspect_ratio: float = Input(1.70, validator = Range(1.5, 2.0))                 # Vertical tail aspect ratio
    ver_tail_volume: float = Input(0.065, validator = Range(0.045, 0.075))                  # Vertical tail volume
    ver_tail_root_chord: float = Input(0.90, validator = GreaterThan(0.0))                  # Vertical tail root chord
    ver_tail_sweep: float = Input(15.0, validator=Range(0.0, 20.0))                         # Leading edge sweep angle [deg]
    ver_tail_avl_n_chordwise = Input(8, validator = GreaterThan(0))                         # Vertical tail chordwise elements for avl
    ver_tail_avl_n_spanwise = Input(12, validator = GreaterThan(0))                         # Vertical tail spanwise elements for avl

    mesh_deflection: float = Input(1e-4, validator = GreaterThan(0.0))

    # Fuselage input parameters
    @Input(validator = Range(0.3, 1.8))
    def fuselage_max_diameter(self) -> float:
        val = 0.0
        if self.occupants == 1:
            val = 0.7
        if self.occupants == 2:
            val = 1.0
        return val
    
    @Input(validator = Range(4.0, 10.0))
    def fuselage_length(self) -> float:
        val = 0.0
        if self.occupants == 1:
            val = 6.5
        if self.occupants == 2:
            val = 8.5
        return val
    
    # Top-level attributes
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
    def cl_cruise_light(self):
        empty_plus_pilot = self.glider_empty_mass + self.current_pilot_mass
        return (empty_plus_pilot*G0)/(0.5*RHO0*self.cruise_speed**2*self.wing_surface_area)
    
    @Attribute
    def cl_cruise_heavy(self):
        return (self.max_to_mass*G0)/(0.5*RHO0*self.cruise_speed**2*self.wing_surface_area)
        
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

    # Fuselage attributes
    @Attribute
    def fuselage_max_radius(self):
        return self.fuselage_max_diameter/2

    # Main wing attributes
    @Attribute
    def wing_span(self):
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
    def wing_surface_area(self):
        return self.max_to_mass / self.wing_loading

    @Attribute
    def wing_position(self):
        return translate(self.position,
                        'x', self.wing_pos_long * self.fuselage_length,
                        'z', self.wing_pos_vert * self.fuselage_max_radius)
    
    @Attribute
    def wing_aspect_ratio(self):
        return self.wing_span**2 / self.wing_surface_area


    # Horizontail tail attributes
    @Attribute
    def hor_tail_tip_chord(self):
        return self.hor_tail_root_chord * self.hor_tail_taper
    
    @Attribute
    def hor_tail_position(self):
        return translate(rotate(self.vert_tail.profiles[-1].position, "x", -90, deg = True),
                         "x", -self.hor_tail_overhang)
    
    @Attribute
    def hor_tail_aspect_ratio(self):
        return 2*self.right_hor_tail.wing_half_aspect_ratio
   
    @Attribute
    def hor_tail_length(self):
        return self.right_hor_tail.x_ac - self.right_wing.x_ac
    
    # Vertical tail attributes
    @Attribute
    def ver_tail_position(self):
        return rotate(translate(self.position, 'x', self.fuselage_length),
                      'x', np.radians(90))
    
    @Attribute
    def ver_tail_length(self):
        quarter_chord_tail = 0.25 * self.ver_tail_root_chord
        return (self.ver_tail_position.x + quarter_chord_tail) - self.right_wing.x_ac
    
    @Attribute
    def dCLv_dBeta(self):
        return (2*np.pi*self.ver_tail_aspect_ratio)/(2 + np.sqrt(self.ver_tail_aspect_ratio**2 + 4))
    
    @Attribute
    def ver_tail_area(self):
        return (self.ver_tail_volume * self.wing_surface_area * self.wing_span) / (self.ver_tail_length * self.dCLv_dBeta)
    
    @Attribute
    def ver_tail_taper(self):
        return self.hor_tail_root_chord / self.ver_tail_root_chord
    
    @Attribute
    def ver_tail_height(self):
        return np.sqrt(self.ver_tail_aspect_ratio * self.ver_tail_area)

    # Parts
    @Part
    def right_wing(self):
        return LiftingSurface(
            name = "Main wing",
            root_chord = self.wing_root_chord,
            span = self.wing_span/2,
            num_sections = 1,                               # Default setting is one wing section, user can add more in GUI
            mesh_deflection = self.mesh_deflection,
            af_cst_order = 5,                               # Default setting, adjustable in GUI under right_wing part
            af_num_points = 200,                            # Default setting, adjustable in GUI under right_wing part 
            af_closed_TE = True,                            # Default setting, adjustable in GUI under right_wing part 
            position = self.wing_position,                  
            has_winglet = self.has_winglet,
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
            n_chordwise = self.wing_avl_n_chordwise,
            chord_spacing = avl.Spacing.cosine,
            n_spanwise = self.wing_avl_n_spanwise,
            span_spacing = avl.Spacing.cosine,
            y_duplicate = self.right_wing.position.point[1],
            sections = [profile.avl_section for profile in self.right_wing.profiles]   
        )

    @Part
    def right_hor_tail(self):
        return LiftingSurface(
            name = "Horizontal tail",
            root_af = "NACA 0010",                              # Hard-coded default value to keep glider class as clean as possible, editable in GUI
            tip_af = "NACA 0010",                               # "
            root_chord = self.hor_tail_root_chord,
            taper = self.hor_tail_taper,
            span = self.hor_tail_span/2,
            twist = 0.0,                                        # Hard-coded default value to keep glider class as clean as possible, editable in GUI
            sweep = 10.0,                                        # "
            sweep_loc = 0.0,                                    # "
            dihedral = 0.0,                                     # "
            incidence_angle = 2.0,                              # "    
            num_sections = 1,                                   # "
            mesh_deflection = self.mesh_deflection,                                 
            af_cst_order = 5,                                   # Hard-coded default value to keep glider class as clean as possible, editable in GUI
            af_num_points = 40,                                 # "
            af_closed_TE = True,                                # "
            position = self.hor_tail_position,
            has_winglet = False                                 # Horizontal and vertical tails almost never have winglets
        )

    @Part
    def left_hor_tail(self):
        return MirroredShape(
            shape_in = self.right_hor_tail,
            reference_point = self.hor_tail_position,
            vector1 = self.position.Vz,
            vector2 = self.position.Vx
        )
    
    @Attribute(in_tree = True)
    def horizontal_tail_avl_surface(self):
        return avl.Surface(
            name = "Horizontal tail",
            n_chordwise = self.hor_tail_avl_n_chordwise,
            chord_spacing = avl.Spacing.cosine,
            n_spanwise = self.hor_tail_avl_n_spanwise,
            span_spacing = avl.Spacing.cosine,
            y_duplicate = self.right_hor_tail.position.point[1],
            sections = [profile.avl_section for profile in self.right_hor_tail.profiles]   
        )

    @Part
    def vert_tail(self):
        return LiftingSurface(
            name = "Vertical tail",
            root_af = "NACA 0012",                              # Hard-coded default value to keep glider class as clean as possible, editable in GUI
            tip_af = "NACA 0012",                               # "
            root_chord = self.ver_tail_root_chord,
            taper = self.ver_tail_taper,
            span = self.ver_tail_height,
            twist = 0.0,                                        # Glider vertical tails do not have twist. Editable in GUI but not relevant
            sweep = self.ver_tail_sweep,                    
            sweep_loc = 0.0,                                    # Hard-coded default value to keep glider class as clean as possible, editable in GUI
            dihedral = 0.0,                                     # Glider vertical tails do not have dihedral (rotation handled by ver_tail_position). Editable in GUI but not relevant
            incidence_angle = 0.0,                              # Glider vertical tails do not have an incidence angle. Editable in GUI but not relevant
            num_sections = 1,                                   # Vertical tails rarely have multiple sections
            mesh_deflection = self.mesh_deflection,                             
            af_cst_order = 5,                                   # Hard-coded default value to keep glider class as clean as possible, editable in GUI
            af_num_points = 40,                                 # "
            af_closed_TE = True,                                # "
            position = self.ver_tail_position,
            has_winglet = False                                 # Vertical tails have no winglets
        )
    
    @Attribute(in_tree = True)
    def vertical_tail_avl_surface(self):
        return avl.Surface(
            name = "Vertical tail",
            n_chordwise = self.ver_tail_avl_n_chordwise,
            chord_spacing = avl.Spacing.cosine,
            n_spanwise = self.ver_tail_avl_n_spanwise,
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

    @action(label = "Make scissor plots")
    def scissor_plot(self):
        scissorplot = ScissorPlot(
            SM = self.SM,
            actual_x_cog = self.glider_x_cog[0],
            fwd_limit_x_cog = self.glider_x_cog[1],
            bwd_limit_x_cog = self.glider_x_cog[2],
            x_LEMAC = self.right_wing.x_LEMAC,
            x_ac = 0.25,                                            # W.r.t mean aerodynamic chord
            s = self.wing_surface_area,
            Sh_S = self.Sh_S,
            l_h = self.hor_tail_length,
            mac = self.right_wing.mean_aerodynamic_chord,
            cm_ac = self.cm_ac,
            cl_a_min_h = CL_MAX_ESTIMATE,
            cl_h = -0.35*self.hor_tail_aspect_ratio**(1/3),
            cl_alpha_h = self.dcl_da_tail,
            cl_alpha_a_min_h = self.dcl_da_wing,
            Vh_V = 1
        )
        scissorplot.plot_scissor_plot()
    
    @action(label = "Size tail")
    def size_tail(self):
        CL_h = -0.35*self.hor_tail_aspect_ratio**(1/3)
        self.Sh_S = size_tail(
            self.glider_x_cog[1], 
            self.glider_x_cog[0],
            self.glider_x_cog[2],
            self.right_wing.mean_aerodynamic_chord,
            self.right_wing.x_LEMAC,
            self.x_np,
            self.cm_ac,
            CL_MAX_ESTIMATE,
            self.dcl_da_tail,
            self.dcl_da_wing,
            CL_h,
            1.0,
            self.hor_tail_length
            )
    
    @Part
    def Q3D_params(self):
        return Q3DData(
            velocity = self.cruise_speed,
            altitude = 1.0
        )

    @action(label = "Run Q3D - speed polar")
    def q3d_data(self):
        """All inputs and results from running Q3D (MATLAB)"""
        res_light = MATLAB_Q3D_ENGINE.run_q3d_cst(
            self.right_wing.q3d_planform_geom,
            self.right_wing.q3d_cst_airfoils,
            self.right_wing.q3d_eta_airfoils,
            matlab.double(self.right_wing.incidence_angle),
            self.Q3D_params.mach_number,
            self.Q3D_params.reynolds_number,
            self.Q3D_params.velocity,
            self.cl_cruise_light,
            self.Q3D_params.altitude,
            self.Q3D_params.density
        )
        res_heavy = MATLAB_Q3D_ENGINE.run_q3d_cst(
            self.right_wing.q3d_planform_geom,
            self.right_wing.q3d_cst_airfoils,
            self.right_wing.q3d_eta_airfoils,
            matlab.double(self.right_wing.incidence_angle),
            self.Q3D_params.mach_number,
            self.Q3D_params.reynolds_number,
            self.Q3D_params.velocity,
            self.cl_cruise_heavy,
            self.Q3D_params.altitude,
            self.Q3D_params.density
        )
        self.q3d_res = [res_light, res_heavy]

    @Attribute
    def cruise_L_D_light(self) -> float:
        result = getattr(self, "q3d_res", None)
        if result is not None:
            return result[0]["CLwing"] / result[0]["CDwing"]
        else:
            return "Evaluate aerodynamics to view property"

    @Attribute
    def cruise_L_D_heavy(self) -> float:
        result = getattr(self, "q3d_res", None)
        if result is not None:
            return result[1]["CLwing"] / result[1]["CDwing"]
        else:
            return "Evaluate aerodynamics to view property"
        
    @Attribute(in_tree = True)
    def avl_surfaces(self):
        return [self.main_wing_avl_surface, self.horizontal_tail_avl_surface, self.vertical_tail_avl_surface]
    
    @Part
    def avl_configuration(self):
        """Configurations are made separately for each Mach number that is provided."""
        return avl.Configuration(name='avl analysis',
                                 reference_area=self.wing_surface_area,
                                 reference_span=self.wing_span,
                                 reference_chord=self.right_wing.mean_aerodynamic_chord,
                                 reference_point=Point(self.right_wing.x_ac, 0.0, self.right_wing.position.z),
                                 surfaces=self.avl_surfaces,
                                 mach= 0.0)
    
    @Attribute
    def avl_find_Cm_ac_settings(self):
        return {'alpha': avl.Parameter(name='alpha',
                                         setting='CL',
                                         value= 0.0)}

    @Attribute
    def avl_dcl_da_settings(self):
        alpha_range = np.arange(0.0, 11.0, 1.0)
        return [{'alpha': aoa} for aoa in alpha_range]
    
    @Part
    def avl_Cm_ac_case(self):
        """avl case definition using the avl_settings dictionary defined above"""
        return avl.Case(name='fixed_cl',  # name _must_ correspond to type of case
                        settings=self.avl_find_Cm_ac_settings)
    
    @Part
    def avl_dcl_da_case(self):
        return avl.Case(quantify = len(self.avl_dcl_da_settings), 
                        name = 'fixed_aoa_' + str(self.avl_dcl_da_settings[child.index]['alpha']),
                        settings = self.avl_dcl_da_settings[child.index])
    
    @Part
    def avl_Cm_ac_analysis(self):
        return avl.Interface(configuration=self.avl_configuration,
                             # note: AVL always expects a list of cases!
                             cases=[self.avl_Cm_ac_case])
    
    @Part
    def avl_dcl_da_analysis(self):
        return avl.Interface(configuration = self.avl_configuration,
                             cases = self.avl_dcl_da_case)
    
    @Attribute
    def dcl_da_wing(self):
        alpha_range = np.arange(0.0, 11.0, 1.0)
        cl_wing = np.array([result['SurfaceForces']['Main Wing']['CL'] for case_nr, result in self.avl_dcl_da_analysis.results.items()])
        return np.rad2deg(np.mean(np.diff(cl_wing)/np.diff(alpha_range)))

    @Attribute
    def dcl_da_tail(self): # Check the convention, what area is used to normalise tail CL. Value seems low
        alpha_range = np.arange(0.0, 11.0, 1.0)
        cl_wing = np.array([result['SurfaceForces']['Horizontal tail']['CL'] for case_nr, result in self.avl_dcl_da_analysis.results.items()])
        return (1/self.Sh_S)*np.rad2deg(np.mean(np.diff(cl_wing)/np.diff(alpha_range)))
    
    @Attribute
    def cm_ac(self):
        cm_ac = np.array([result['SurfaceForces']['Main Wing']['Cm'] for case_nr, result in self.avl_Cm_ac_analysis.results.items()])
        return np.mean(cm_ac)
    
    @Attribute
    def x_np(self):
        x_np = np.array([result['StabilityDerivatives']['Xnp'] for case_nr, result in self.avl_dcl_da_analysis.results.items()])
        return np.mean(x_np)
    
    @Part
    def step_writer(self):
        return STEPWriter(trees = [self],
                          filename = Path.cwd() / "output" / "glider.stp")
    
    @action(label = "Write STEP file")
    def write_step_file(self):
        self.step_writer.write()

if __name__ == '__main__':
    from parapy.gui import display
    glider = Glider(name="glider")
    display(glider)