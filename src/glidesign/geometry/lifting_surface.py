# Python native imports
from typing import Optional

# Python third party imports
import numpy as np
import matlab

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, LoftedSolid, Position, Point, Orientation, Vector
from parapy.core import Input, Attribute, Part, Base
from parapy.core.validate import OneOf, Range, GreaterThan, Validator, GreaterThanOrEqualTo, LessThan

# Custom imports
from ..core import airfoil_found
from .airfoil import Airfoil
from .ref_frame import Frame
from ..external import MATLAB_Q3D_ENGINE

class LiftingSection(GeomBase):
    """
    Serves as a parameter container for multi section LiftingSurface
    
    """
    idx: int = Input(validator = GreaterThanOrEqualTo(0))                   # Section index, later used in LiftingSurface
    root_af: str = Input(validator = airfoil_found)                         # ID of the section root airfoil
    tip_af: str = Input(validator = airfoil_found)                          # ID of the section tip airfoil
    root_chord: float = Input(0.5, validator = GreaterThan(0.0))            # Root chord of the section [m]
    tip_chord: float = Input(0.2, validator = GreaterThan(0.0))             # Tip chord of the section [m]              
    span: float = Input(7.3, validator = GreaterThan(0.0))                  # Section span [m]
    twist: float = Input(0.0)                                               # Section twist angle [deg]
    dihedral: float = Input(0.0)                                            # Section dihedral angle [deg]
    sweep: float = Input(0.0)                                               # Section sweep angle [deg]
    sweep_loc: float = Input(0.25, validator = Range(0.0, 1.0))  # Chord normalised coordinate around which sweep is applied

    @Validator
    def validate_tip_chord(self):
        if self.tip_chord >= self.root_chord:
            raise ValueError(f"Tip chord must be smaller than root chord")
        return True

    @Attribute
    def taper_ratio(self):
        return self.tip_chord/self.root_chord
    
    @Attribute
    def area(self):
        return 0.5 * (self.root_chord + self.tip_chord) * self.span
    
    @Attribute
    def section_mean_aerodynamic_chord(self):
        return (2/3)*self.root_chord*((1+self.taper_ratio+self.taper_ratio**2)/(1+self.taper_ratio))
    
class Q3DData(Base):

    mach = Input(0.0, validator = GreaterThanOrEqualTo(0.0))            # Mach number for adding compressibility effects to Q3D [-]
    reynolds = Input(30.0, validator = GreaterThan(0.0))                # Reynolds number for viscous analysis of Q3D [-]
    velocity = Input(30.0, validator = GreaterThan(0.0))                # Velocity at which Q3D analysis is run [m/s]
    alpha = Input(2.0, validator = LessThan(14.0))                      # Angle of attack at which Q3D is run [deg]
    altitude = Input(500.0, validator = GreaterThan(0.0))               # Altitude for Q3D [m]  
    density = Input(1.225, validator = GreaterThan(0.0))                # Air density for Q3D analysis [kg/m^3]

class LiftingSurface(LoftedSolid):

    name: str = Input()                                                 # Name of the part (e.g. Wing, Horizontal Tail)
    incidence_angle: float = Input(0.0)                                 # Incidence angle of the lifting surface [deg]
    sections: list[LiftingSection] = Input([                            # List of sections, can be edited on runtime
        LiftingSection(
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
            )
        ])
    mesh_deflection: float = Input(1e-4)                                # Parameter for LoftedSolid superclass
    af_cst_order: int = Input(5)                                        # Polynomial order to use for CST representation of airfoils
    af_num_points: int = Input(30)                                      # Number of points to use for NACA airfoil generation
    af_closed_TE: bool = Input(True)                                    # Closed trailing edges for surface airfoils
    Q3D_params: Q3DData = Input(Q3DData())

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
    def num_sections(self) -> int:
        return len(self.sections)
    
    @Attribute
    def profiles(self):
        profile_list = [Airfoil(
            airfoil_name = self.sections[0].root_af,
            chord = self.sections[0].root_chord,
            cst_poly_order = self.af_cst_order,
            num_points = self.af_num_points,
            position = rotate(self.position, 'y', self.incidence_angle, deg=True)
        )]
        for i in range(self.num_sections):
            next_pos = profile_list[-1].position
            next_pos = rotate(next_pos, "y", self.sections[i].twist, deg = True)
            next_pos = translate(next_pos,
                                "x", self.sections[i].span * np.tan(np.deg2rad(self.sections[i].sweep)) + 
                                (self.sections[i].root_chord - self.sections[i].tip_chord) * self.sections[i].sweep_loc,
                                "y", self.sections[i].span,
                                "z", self.sections[i].span * np.tan(np.deg2rad(self.sections[i].dihedral))
                            )
            profile_list.append(Airfoil(
                airfoil_name = self.sections[i].tip_af,
                chord = self.sections[i].tip_chord,
                cst_poly_order = self.af_cst_order,
                num_points = self.af_num_points,
                position = next_pos
            ))
        return profile_list

#Sorry als ik dit een beetje verpest heb, maar heb voor mezelf ff hier de full parameters toegevoegd (had ze nodig voor scissor_plot.py

    @Attribute
    def semi_span(self):
        return sum([sec.span for sec in self.sections])

    @Attribute
    def span(self):
        return self.semi_span * 2
    
    @Attribute
    def half_area(self):
        return sum([sec.area for sec in self.sections])

    @Attribute
    def full_area(self):
        return self.half_area * 2
    
    @Attribute
    def wing_half_aspect_ratio(self):
        return self.semi_span*self.semi_span / self.half_area

    @Attribute
    def wing_full_aspect_ratio(self):
        return self.wing_half_aspect_ratio * 2
    
    @Attribute
    def mean_aerodynamic_chord(self):
        return sum([sec.section_mean_aerodynamic_chord*sec.span for sec in self.sections])/self.semi_span

    a = wing_half_aspect_ratio

    @Attribute
    def q3d_planform_geom(self) -> matlab.double:
        root_x, root_y, root_z = self.profiles[0].position.location
        planform_geom = [[root_x, root_y, root_z, self.profiles[0].chord, 0.0]]
        for i in range(1, self.num_sections):
            tip_x, tip_y, tip_z = self.profiles[i].position.location
            twist_angle = self.profiles[i].position.orientation
            planform_geom.append([tip_x, tip_y, tip_z, self.profiles[i].chord, twist_angle])
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
        eta = [0.0]
        for i in range(1, self.num_sections + 1):
            eta.append(self.profiles[i].position.location.y / self.semi_span)
        return matlab.double(eta)
    
    @Attribute
    def q3d_data(self):
        """All inputs and results from running Q3D (MATLAB)"""
        return MATLAB_Q3D_ENGINE.run_q3d_cst(
            self.q3d_planform_geom,
            self.q3d_cst_airfoils,
            self.q3d_eta_airfoils,
            matlab.double(self.incidence_angle),
            self.Q3D_params.mach,
            self.Q3D_params.reynolds,
            self.Q3D_params.velocity,
            self.Q3D_params.alpha,
            self.Q3D_params.altitude,
            self.Q3D_params.density
        )
    
    @Attribute
    def q3d_res(self) -> dict:
        """q3d results"""
        return self.q3d_data[0]

    @Attribute
    def q3d_ac(self) -> dict:
        """q3d inputs"""
        return self.q3d_data[1]

    @Attribute
    def wing_cl(self) -> float:
        return self.q3d_res["CLwing"]

    @Attribute
    def wing_cd(self) -> float:
        return self.q3d_res["CDwing"]
    
    @Part
    def frame(self):
        """to visualize the given lifting surface reference frame"""
        return Frame(pos=self.position,
                     hidden=False)