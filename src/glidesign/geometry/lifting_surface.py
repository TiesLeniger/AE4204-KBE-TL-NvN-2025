# Python native imports
from typing import Optional

# Python third party imports
import numpy as np
import matlab

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, LoftedSolid, Position
from parapy.core import Input, Attribute, Part
from parapy.core.validate import OneOf, Range, GreaterThan, Validator, GreaterThanOrEqualTo

# Custom imports
from ..core import airfoil_found, validate_equal_length_lists
from .airfoil import Airfoil
from .ref_frame import Frame

class LiftingSection(GeomBase):
    
    section_idx: int = Input(0, validator = GreaterThanOrEqualTo(0))
    previous_section: Optional["LiftingSection"] = Input(None)          # Previous section (used for positioning)
    parent_surface: "LiftingSurface" = Input()

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
    
    mesh_deflection: float = Input(1e-4)

    @Attribute
    def root_position(self):
        if self.section_idx == 0:
            return self.parent_surface.position
        elif self.previous_section:
            return self.previous_section.tip_position
        else:
            raise ValueError(f"Section index is larger than 0 ({self.section_idx}) but a previous section was not provided")
    
    @Attribute
    def tip_position(self):
        rotated_pos = rotate(self.root_position, "y", self.twist, deg = True)
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

class LiftingSurface(GeomBase):
    # Global parameters
    name: str = Input()
    position: Position = Input(Position(0.0, 0.0, 0.0))
    orientation: tuple[float, float, float] = Input((0.0, 0.0, 0.0))
    surface_root_chord: float = Input(0.5, validator = GreaterThan(0.0))

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
    
    @Attribute
    def num_sections(self) -> int:
        return len(self.sec_spans)

    @Attribute
    def sec_root_chords(self) -> list[float]:
        sec_root_chords = [self.surface_root_chord]
        for taper in self.sec_tapers[:-1]:
            sec_root_chords.append(sec_root_chords[-1] * taper)
        return sec_root_chords

    @Attribute
    def sec_tip_chords(self) -> list[float]:
        return [root * taper for root, taper in zip(self.sec_root_chords, self.sec_tapers)]
    
    @Part
    def sections(self):
        def builder(i):
            return LiftingSection(
                section_idx = i,
                previous_section = self.sections[i-1] if i>0 else None,
                parent_surface = self,
                root_airfoil_id = self.root_airfoil_ids[i],
                tip_airfoil_id = self.tip_airfoil_ids[i],
                root_chords = self.sec_root_chords[i],
                taper_ratio = self.sec_tapers[i],
                span = self.sec_spans[i],
                twist = self.sec_twists[i],
                dihedral=self.sec_dihedrals[i],
                sweep=self.sec_sweeps[i],
                sweep_loc=self.sec_sweeplocs[i],
                control_surface="None",                 # TODO: add controlsurface functionality
                mesh_deflection=self.mesh_deflection 
            )
        return builder, self.num_sections
    
    @Part
    def lofted_surface(self):
        return LoftedSolid(
            profiles=[s.root_airfoil for s in self.sections] + [self.sections[-1].tip_airfoil],
            mesh_deflection=self.mesh_deflection
        )

if __name__ == '__main__':
    from parapy.gui import display
    wing = LiftingSurface(name="wing_test")
    display(wing)