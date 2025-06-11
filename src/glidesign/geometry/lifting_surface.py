# Python native imports
from typing import Optional

# Python third party imports
import numpy as np
import matlab

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, LoftedSolid, Position, Point, Orientation, Vector
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range, GreaterThan, Validator, GreaterThanOrEqualTo

# Custom imports
from ..core import airfoil_found
from .airfoil import Airfoil
from .ref_frame import Frame

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
    sweep_loc: float = Input(0.25, validator = Range(0.0, 1.0))             # Chord normalised coordinate around which sweep is applied

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

    @Attribute
    def semi_span(self):
        return sum([sec.span for sec in self.sections])
    
    @Attribute
    def area(self):
        return sum([sec.area for sec in self.sections])
    
    @Attribute
    def wing_half_aspect_ratio(self):
        return self.semi_span*self.semi_span / self.area
    
    @Attribute
    def mean_aerodynamic_chord(self):
        return sum([sec.section_mean_aerodynamic_chord*sec.span for sec in self.sections])/self.semi_span