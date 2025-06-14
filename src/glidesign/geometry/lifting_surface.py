# Python native imports
from typing import Optional
from copy import deepcopy

# Python third party imports
import numpy as np
import matlab

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, LoftedSolid, Position, Point, Orientation, Vector
from parapy.core import Input, Attribute, Part, Base, action, PartList
from parapy.core.validate import OneOf, Range, GreaterThan, Validator, GreaterThanOrEqualTo, LessThan

# Custom imports
from ..core import airfoil_found, convert_matlab_dict
from .airfoil import Airfoil
from .ref_frame import Frame
from ..external import MATLAB_Q3D_ENGINE, Q3DData

class LiftingSection(GeomBase):
    """
    Serves as a parameter container for multi section LiftingSurface
    """
    idx: int = Input(validator = GreaterThanOrEqualTo(0))                   # Section index, later used in LiftingSurface
    root_af_id: str = Input(validator = airfoil_found)                      # ID of the section root airfoil
    tip_af_id: str = Input(validator = airfoil_found)                       # ID of the section tip airfoil
    root_chord: float = Input(0.5, validator = GreaterThan(0.0))            # Root chord of the section [m]
    tip_chord: float = Input(0.2, validator = GreaterThan(0.0))             # Tip chord of the section [m]              
    span: float = Input(7.3, validator = GreaterThan(0.0))                  # Section span [m]
    twist: float = Input(0.0)                                               # Section twist angle [deg]
    dihedral: float = Input(0.0)                                            # Section dihedral angle [deg]
    sweep: float = Input(0.0)                                               # Section sweep angle [deg]
    sweep_loc: float = Input(0.25, validator = Range(0.0, 1.0))             # Chord normalised coordinate around which sweep is applied
    af_cst_order: int = Input(4, validator = GreaterThan(0))                # CST order for the airfoils
    af_num_points: int = Input(200, validator = GreaterThan(0))             # Number of points on the airfoils (for NACA generation)
    mesh_deflection: float = Input(1e-4, validator = GreaterThan(0.0))      # mesh deflection parameter for lofted surface

    @Validator
    def validate_tip_chord(self):
        if self.tip_chord >= self.root_chord:
            raise ValueError(f"Tip chord must be smaller than root chord")
        return True
    
    @Part
    def root_airfoil(self):
        return Airfoil(
            airfoil_name = self.root_af_id,
            chord = self.root_chord,
            cst_poly_order = self.af_cstorder,
            num_points = self.af_numpoints,
            mesh_deflection = self.mesh_deflection,
            position = self.position
        )
    
    @Attribute
    def tip_position(self):
        tip_position = rotate(self.position, "y", self.twist, deg = True)
        tip_position = translate(tip_position,
                            "x", self.span * np.tan(np.deg2rad(self.sweep)) + (self.root_chord - self.tip_chord) * self.sweep_loc,
                            "y", self.span,
                            "z", self.span * np.tan(np.deg2rad(self.dihedral))
                        )
        return tip_position
    
    @Part
    def tip_airfoil(self):
        return Airfoil(
            airfoil_name = self.tip_af_id,
            chord = self.tip_chord,
            cst_poly_order = self.af_cstorder,
            num_points = self.af_numpoints,
            mesh_deflection = self.mesh_deflection,
            position = self.tip_position
        )
        
    @Attribute
    def taper_ratio(self):
        return self.tip_chord / self.root_chord
    
    @Attribute
    def area(self):
        return 0.5 * (self.root_chord + self.tip_chord) * self.span
    
    @Attribute
    def section_mean_aerodynamic_chord(self):
        return (2/3)*self.root_chord*((1+self.taper_ratio+self.taper_ratio**2)/(1+self.taper_ratio))
    
class LiftingSurface(LoftedSolid):

    name: str = Input()                                                 # Name of the part (e.g. Wing, Horizontal Tail)
    root_af: str = Input(validator = airfoil_found)
    tip_af: str = Input(validator = airfoil_found)
    root_chord: float = Input(0.5, validator = GreaterThan(0.0))
    taper: float = Input(0.5, validator = Range(0.0, 1.0))
    span: float = Input(7.5, validator = GreaterThan(0.0))
    twist: float = Input(0.0)
    sweep: float = Input(0.0)
    sweep_loc: float = Input(0.25, validator = Range(0.0, 1.0))
    dihedral: float = Input(0.0)
    incidence_angle: float = Input(0.0)                                 # Incidence angle of the lifting surface [deg]
    mesh_deflection: float = Input(1e-4)                                # Parameter for LoftedSolid superclass
    af_cst_order: int = Input(4)                                        # Polynomial order to use for CST representation of airfoils
    af_num_points: int = Input(30)                                      # Number of points to use for NACA airfoil generation
    af_closed_TE: bool = Input(True)                                    # Closed trailing edges for surface airfoils
    ruled: bool = Input(True)

    has_winglet: bool = Input(False)                                                    # Boolean for adding a winglet
    winglet_length: float = Input(0.0, validator = Range(0.0, 1.0, incl_min = False))   # Winglet length in [m]
    winglet_cant: float = Input(5.0, validator = Range(0.0, 90.0))                      # Cant angle of the winglet [deg] (90 deg means wing extension)
    winglet_toe: float = Input(1.0, Range(-5.0, 5.0))                                   # Toe angle of the winglet [deg]
    winglet_sweep: float = Input(30.0, Range(0.0, 50.0))                                # Leading edge sweep of the winglet [deg]
    winglet_tip_af: str = Input('NACA0010', validator = airfoil_found)                  # Airfoil profile of the winglet
    winglet_taper: float = Input(0.5, Range(0.1, 1.0))                                  # Winglet taper ratio
    
    @Input(validator = Range(1, 10))
    def num_sections(self):
        return len(self.sections)

    @Validator
    def validate_section_airfoils(self):
        if len(self.sections) > 1:
            for i in range(self.num_sections - 1):
                if self.sections[i].tip_af != self.sections[i+1].root_af:
                    raise ValueError(
                        f"The tip airfoil of section {i} must be equal to the root airfoil of section {i+1}"
                        )
        return True
    
    @Attribute
    def sections(self):
        sections = [LiftingSection(
            idx = i,
            root_af_id = self.root_af,
            tip_af_id = self.tip_af,
            root_chord = self.root_chord,
            tip_chord = self.root_chord * self.taper**(1/self.num_sections),
            span = self.span / self.num_sections,
            twist = self.twist / self.num_sections,
            dihedral = self.dihedral,
            sweep = self.sweep,
            sweep_loc = self.sweep_loc,
            af_cst_order = self.af_cst_order,
            af_num_points = self.af_num_points,
            mesh_deflection = self.mesh_deflection,
            position = self.position if i == 0 else sections[i-1].tip_position
            ) for i in range(self.num_sections)]
        return sections

    @Part
    def profiles(self):
        profiles = [self.sections[0].root_airfoil]
        for i in range(self.num_sections):
            profiles.append(self.sections[i].tip_airfoil)

    @Part
    def lofted_solid(self):
        return LoftedSolid(profiles = self.profiles, hidden = False)
    
    @Part
    def winglet(self):
        if self.has_winglet:
            winglet_pos = self.profiles[-1].position
            winglet_pos = rotate(winglet_pos, 
                                 "x", 90 - self.winglet_cant,
                                 "z", self.winglet_toe)
            winglet = LiftingSection(
                idx = self.num_sections,
                root_af_id = self.sections[-1].tip_af_id,
                root_chord = self.sections[-1].tip_chord,
                span = self.winglet_length,
                twist = 0.0,
                dihedral = 0.0,
                sweep = self.winglet_sweep,
                sweep_loc = 0.0,
                af_cst_order = self.af_cst_order,
                af_num_points = self.af_num_points,
                mesh_deflection = self.mesh_deflection,
                position = winglet_pos
            )
            self.sections.append(winglet)
            return winglet
        else:
            return None

    @Attribute
    def semi_span(self):
        return sum([sec.span for sec in self.sections])
    
    @Attribute
    def half_area(self):
        return sum([sec.area for sec in self.sections])
    
    @Attribute
    def wing_half_aspect_ratio(self):
        return self.semi_span*self.semi_span / self.half_area
    
    @Attribute
    def mean_aerodynamic_chord(self):
        return sum([sec.section_mean_aerodynamic_chord*sec.span for sec in self.sections])/self.semi_span

    @Attribute
    def q3d_planform_geom(self) -> matlab.double:
        planform_geom = []
        for i in range(len(self.profiles)):
            x, y, z = self.profiles[i].position.location
            twist = 0.0 if i == 0 else self.sections[i-1].twist
            planform_geom.append([x, y, z, self.profiles[i].chord, twist])
        return matlab.double(planform_geom)
    
    @Attribute
    def q3d_cst_airfoils(self) -> matlab.double:
        cst_coeffs = []
        cst_coeffs.append(self.profiles[0].cst_coeff_u + self.profiles[0].cst_coeff_l)
        for i in range(1, self.num_sections + 1):
            cst_coeffs.append(self.profiles[i].cst_coeff_u + self.profiles[i].cst_coeff_l)
        cst_coeffs = matlab.double(cst_coeffs)

        return cst_coeffs
    
    @Attribute
    def q3d_eta_airfoils(self) -> matlab.double:
        eta = [[0.0]]
        for i in range(1, self.num_sections + 1):
            eta.append([self.profiles[i].position.location.y / self.semi_span])
        return matlab.double(eta)

    @Part
    def Q3D_params(self):
        return Q3DData()
    
    @action(label = "Run Q3D")
    def q3d_data(self):
        """All inputs and results from running Q3D (MATLAB)"""
        self.q3d_res = MATLAB_Q3D_ENGINE.run_q3d_cst(
            self.q3d_planform_geom,
            self.q3d_cst_airfoils,
            self.q3d_eta_airfoils,
            matlab.double(self.incidence_angle),
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
    
    @Part
    def frame(self):
        """to visualize the given lifting surface reference frame"""
        return Frame(pos=self.position,
                     hidden=False)