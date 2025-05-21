# Python native imports
from typing import Optional

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, MirroredShape, LoftedSolid, Position
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range, GreaterThan, Validator

# Custom imports
from ..core import airfoil_found, validate_equal_length_lists
from .airfoil import Airfoil
from .ref_frame import Frame

class LiftingSection(LoftedSolid):

    position: Optional[Position] = Input(Position(0.0, 0.0, 0.0))
    orientation: Optional[tuple[float, float, float]] = Input((0.0, 0.0, 0.0))

    root_airfoil_id: str = Input(validator = airfoil_found)             # Inner airfoil name
    tip_airfoil_id: str = Input(validator = airfoil_found)              # Outer airfoil name

    root_chord: float = Input(validator = GreaterThan(0))               # Root chord [m]
    taper_ratio: float = Input(validator = Range(0.0, 1.0))             # Taper ratio of the section
    span: float = Input(7.5, validator = GreaterThan(0))                # Section span [m]

    twist: float = Input(0.0, validator = Range(-5.0, 5.0))             # Section twist (tip w.r.t root) [deg]
    dihedral: float = Input(0.0, validator = Range(-3.0, 5.0))          # Section dihedral [deg]
    sweep: float = Input(0.0, validator = Range(-10.0, 20.0))           # Section sweep [deg]
    sweep_loc: float = Input(0.25, validator = Range(0.0, 1.0))         # Chord normalized sweep location of the section (0 is LE, 1 is TE)
    
    control_surface: str = Input("None", OneOf(                         # Control surface on the section
        "None", "Aileron", "Flaperon", "Flap", "Airbrake"))
    
    previous_section: Optional["LiftingSection"] = Input(None)          # Previous section (used for positioning)
    
    mesh_deflection: float = Input(1e-4)
    
    @Attribute
    def profiles(self):
        return [self.root_airfoil, self.tip_airfoil]

    @Attribute
    def root_position(self):
        if self.previous_section is None:
            return rotate(self.position,
                          "x", self.orientation[0],
                          "y", self.orientation[1],
                          "z", self.orientation[2],
                          deg = True
                        )
        else:
            return self.previous_section.tip_position
    
    @Attribute
    def tip_position(self):
        rotated_pos = rotate(self.root_airfoil.position, "y", self.twist, deg = True)
        sweep_x = self.span * np.tan(np.deg2rad(self.sweep)) + (self.root_chord - self.tip_chord) * self.sweep_loc
        translated_pos = translate(rotated_pos,
                                   "x", sweep_x,
                                   "y", self.span,
                                   "z", self.span * np.tan(np.deg2rad(self.dihedral)))
        
        return translated_pos

    @Attribute
    def tip_chord(self):
        return self.root_chord * self.taper_ratio

    @Part
    def root_airfoil(self):
        return Airfoil(airfoil_name = self.root_airfoil_id,
                       chord = self.root_chord,
                       position = self.root_position,
                       )
    
    @Part
    def tip_airfoil(self):
        return Airfoil(airfoil_name = self.tip_airfoil_id,
                       chord = self.tip_chord,
                       position = self.tip_position,
                       )

class LiftingSurface(LoftedSolid):
    # Global parameters
    name: str = Input()
    position: Position = Input(Position(0.0, 0.0, 0.0))
    orientation: tuple[float, float, float] = Input((0.0, 0.0, 0.0))

    # Section inputs
    root_airfoil_ids: list[str] = Input(["nlf1-0015.dat"])
    tip_airfoil_ids: list[str] = Input(["nlf1-0015.dat"])
    sec_spans: list[float] = Input([7.3])                               # [m]
    sec_twists: list[float] = Input([0.0])                              # [deg]    
    sec_dihedrals: list[float] = Input([0.0])                           # [deg]
    sec_sweeps: list[float] = Input([0.0])                              # [deg]
    sec_sweeplocs: list[float] = Input([0.25])                          # [-]
    sec_tapers: list[float] = Input([0.7])                              # [-]
    
    mesh_deflection: float = Input(1e-4)

    @Validator
    def validate_section_input_lengths(self):
        validate_equal_length_lists(self, [
            "root_airfoil_ids",
            "tip_airfoil_ids",
            "sec_spans",
            "sec_twists",
            "sec_dihedrals",
            "sec_sweeps",
            "sec_sweeplocs",
            "sec_tapers"
        ])

    @Part
    def sections(self):
        previous = None
        sections = []
        n = len(self.root_airfoil_ids)
        for i in range(n):
            if i == 0:
                root_position = self.position
            else:
                root_position = None
            section = LiftingSection(
                
            )