# Python native imports
import numpy as np
import matlab

# ParaPy imports
from parapy.geom import GeomBase, translate, rotate, LoftedSolid
from parapy.core import Input, Attribute, Part, action, child
from parapy.core.validate import Range, GreaterThan, Validator

# Custom imports
from ..core import airfoil_found
from .airfoil import Airfoil
from .ref_frame import Frame

class LiftingSection(GeomBase):
    """
    A class that represents a lifting section of a lifting surface (e.g., a wing). The class 
    defines the geometric properties of the section, including the root and tip airfoils, chord 
    lengths, span, twist, dihedral, and sweep, and computes important aerodynamic properties.

    Attributes:
        root_af_id (str): The ID of the root airfoil.
        tip_af_id (str): The ID of the tip airfoil.
        root_chord (float): The root chord length of the section.
        tip_chord (float): The tip chord length of the section.
        span (float): The span of the section.
        twist (float): The twist angle of the section.
        dihedral (float): The dihedral angle of the section.
        sweep (float): The sweep angle of the section.
        sweep_loc (float): The normalized chord location at which the sweep is applied.
        af_cst_order (int): The CST order for the airfoils.
        af_num_points (int): The number of points for airfoil generation.
        mesh_deflection (float): The mesh deflection parameter for the lofted surface.
    """

    root_af_id: str = Input(validator=airfoil_found)  # Root airfoil ID
    tip_af_id: str = Input(validator=airfoil_found)   # Tip airfoil ID
    root_chord: float = Input(0.9, validator=GreaterThan(0.0))  # Root chord length
    tip_chord: float = Input(0.5, validator=GreaterThan(0.0))   # Tip chord length
    span: float = Input(7.5, validator=GreaterThan(0.0))        # Span of the section
    twist: float = Input(0.0)  # Twist angle of the section [deg]
    dihedral: float = Input(0.0)  # Dihedral angle of the section [deg]
    sweep: float = Input(0.0)  # Sweep angle of the section [deg]
    sweep_loc: float = Input(0.25, validator=Range(0.0, 1.0))  # Sweep location
    af_cst_order: int = Input(4, validator=GreaterThan(0))  # CST order for airfoils
    af_num_points: int = Input(200, validator=GreaterThan(0))  # Number of airfoil points
    mesh_deflection: float = Input(1e-4, validator=GreaterThan(0.0))  # Mesh deflection parameter

    @Validator
    def validate_tip_chord(self):
        """
        Validates that the tip chord is smaller than the root chord.

        Returns:
            bool: True if the tip chord is smaller than the root chord, otherwise raises a ValueError.
        """
        if self.tip_chord >= self.root_chord:
            raise ValueError(f"Tip chord must be smaller than root chord")
        return True
    
    @Part
    def root_airfoil(self):
        """
        Returns the root airfoil object for this section, based on the provided root airfoil ID.
        
        Returns:
            Airfoil: The root airfoil of the section.
        """
        return Airfoil(
            airfoil_name=self.root_af_id,
            chord=self.root_chord,
            cst_poly_order=self.af_cst_order,
            num_points=self.af_num_points,
            mesh_deflection=self.mesh_deflection,
            position=self.position
        )
    
    @Attribute
    def tip_position(self):
        """
        Calculates and returns the position of the tip of the section, considering twist, dihedral, and sweep.

        Returns:
            Point: The position of the tip of the lifting section.
        """
        tip_position = rotate(self.position, "y", self.twist, deg=True)
        tip_position = translate(tip_position,
                                 "x", self.span * np.tan(np.deg2rad(self.sweep)) + (self.root_chord - self.tip_chord) * self.sweep_loc,
                                 "y", self.span,
                                 "z", self.span * np.tan(np.deg2rad(self.dihedral)))
        return tip_position
    
    @Part
    def tip_airfoil(self):
        """
        Returns the tip airfoil object for this section, based on the provided tip airfoil ID.
        
        Returns:
            Airfoil: The tip airfoil of the section.
        """
        return Airfoil(
            airfoil_name=self.tip_af_id,
            chord=self.tip_chord,
            cst_poly_order=self.af_cst_order,
            num_points=self.af_num_points,
            mesh_deflection=self.mesh_deflection,
            position=self.tip_position
        )
        
    @Attribute
    def taper_ratio(self):
        """
        Computes and returns the taper ratio of the section (tip chord / root chord).
        
        Returns:
            float: The taper ratio of the section.
        """
        return self.tip_chord / self.root_chord
    
    @Attribute
    def area(self):
        """
        Computes and returns the area of the section using the formula for the area of a trapezoid.
        
        Returns:
            float: The area of the section.
        """
        return 0.5 * (self.root_chord + self.tip_chord) * self.span
    
    @Attribute
    def section_mean_aerodynamic_chord(self):
        """
        Computes and returns the mean aerodynamic chord (MAC) of the section.
        
        Returns:
            float: The mean aerodynamic chord of the section.
        """
        return (2 / 3) * self.root_chord * ((1 + self.taper_ratio + self.taper_ratio**2) / (1 + self.taper_ratio))
    
class LiftingSurface(LoftedSolid):
    """
    A class representing a lifting surface, such as a wing. The class defines the geometric properties 
    of the surface, such as the root and tip airfoils, chord lengths, span, twist, and dihedral, and 
    generates a 3D solid representation of the surface.

    Attributes:
        name (str): The name of the part (e.g. Wing, Horizontal Tail).
        root_af (str): The ID of the root airfoil.
        tip_af (str): The ID of the tip airfoil.
        root_chord (float): The root chord length of the surface.
        taper (float): The taper ratio of the surface.
        span (float): The span of the surface.
        twist (float): The twist angle of the surface.
        sweep (float): The sweep angle of the surface.
        sweep_loc (float): The normalized chord location at which the sweep is applied.
        dihedral (float): The dihedral angle of the surface.
        incidence_angle (float): The incidence angle of the lifting surface.
        mesh_deflection (float): The mesh deflection parameter.
        af_cst_order (int): The CST order for the airfoils.
        af_num_points (int): The number of points for airfoil generation.
        af_closed_TE (bool): Whether to use closed trailing edges for the airfoils.
        ruled (bool): Whether the surface is ruled.
        has_winglet (bool): Whether the surface has a winglet.
    """
    
    # Input design parameters
    name: str = Input()  # Name of the part (e.g. Wing, Horizontal Tail)
    root_af: str = Input('NACA 65210', validator=airfoil_found)
    tip_af: str = Input('NACA 65210', validator=airfoil_found)
    root_chord: float = Input(0.9, validator=GreaterThan(0.0))
    taper: float = Input(0.4, validator=Range(0.0, 1.0))
    span: float = Input(7.5, validator=GreaterThan(0.0))
    twist: float = Input(0.0)
    sweep: float = Input(2.0)
    sweep_loc: float = Input(0.25, validator=Range(0.0, 1.0))
    dihedral: float = Input(2.0)
    incidence_angle: float = Input(0.0)  # Incidence angle of the lifting surface [deg]
    mesh_deflection: float = Input(1e-4)  # Parameter for LoftedSolid superclass
    af_cst_order: int = Input(4)  # Polynomial order to use for CST representation of airfoils
    af_num_points: int = Input(30)  # Number of points to use for NACA airfoil generation
    af_closed_TE: bool = Input(True)  # Closed trailing edges for surface airfoils
    ruled: bool = Input(True)

    # Winglet-related parameters
    has_winglet: bool = Input(False)
    winglet_length: float = Input(0.3, validator=Range(0.0, 1.0, incl_min=False))  # Winglet length in [m]
    winglet_cant: float = Input(5.0, validator=Range(0.0, 90.0))  # Cant angle of the winglet [deg]
    winglet_toe: float = Input(1.0, Range(-5.0, 5.0))  # Toe angle of the winglet [deg]
    winglet_sweep: float = Input(30.0, Range(0.0, 50.0))  # Leading edge sweep of the winglet [deg]
    winglet_tip_af: str = Input('NACA0010', validator=airfoil_found)  # Airfoil profile of the winglet
    winglet_taper: float = Input(0.5, Range(0.1, 1.0))  # Winglet taper ratio
    
    @Input(validator=Range(1, 10))
    def num_sections(self):
        return len(self.sections)

    @Validator
    def validate_section_airfoils(self):
        """
        Validates that the tip airfoil of each section matches the root airfoil of the next section.
        
        Returns:
            bool: True if the validation passes, otherwise raises a ValueError.
        """
        if len(self.sections) > 1:
            for i in range(self.num_sections - 1):
                if self.sections[i].tip_af != self.sections[i+1].root_af:
                    raise ValueError(
                        f"The tip airfoil of section {i} must be equal to the root airfoil of section {i+1}"
                    )
        return True
    
    @Part
    def sections(self):
        """
        Generates the sections of the lifting surface based on the input parameters.

        Returns:
            LiftingSection: A list of sections representing the lifting surface.
        """
        return LiftingSection(
            quantify=self.num_sections,
            root_af_id=self.root_af if child.index == 0 else child.previous.tip_af_id,
            tip_af_id=self.tip_af,
            root_chord=self.root_chord,
            tip_chord=self.root_chord * self.taper**(1/self.num_sections),
            span=self.span / self.num_sections,
            twist=self.twist / self.num_sections,
            dihedral=self.dihedral,
            sweep=self.sweep,
            sweep_loc=self.sweep_loc,
            af_cst_order=self.af_cst_order,
            af_num_points=self.af_num_points,
            mesh_deflection=self.mesh_deflection,
            position=self.position if child.index == 0 else child.previous.tip_position
        )
    
    @Attribute
    def winglet_pos(self):
        """
        Calculates the position of the winglet based on the input cant and toe angles.
        
        Returns:
            Point: The position of the winglet.
        """
        winglet_pos = rotate(self.sections[-1].tip_position,
                             "x", 90 - self.winglet_cant,
                             deg=True)
        winglet_pos = rotate(winglet_pos,
                             "z", self.winglet_toe,
                             deg=True)
        return winglet_pos
    
    @Part
    def winglet(self):
        """
        Generates the winglet as a part of the lifting surface.

        Returns:
            LiftingSection: The winglet section.
        """
        return LiftingSection(
            suppress=not self.has_winglet,
            root_af_id=self.sections[-1].tip_af_id,
            tip_af_id=self.winglet_tip_af,
            root_chord=self.sections[-1].tip_chord,
            tip_chord=self.winglet_taper * self.sections[-1].tip_chord,
            span=self.winglet_length,
            twist=0.0,
            dihedral=0.0,
            sweep=self.winglet_sweep,
            sweep_loc=0.0,
            af_cst_order=self.af_cst_order,
            af_num_points=self.af_num_points,
            mesh_deflection=self.mesh_deflection,
            position=self.winglet_pos
        )

    @Attribute
    def profiles(self):
        """
        Returns a list of airfoil profiles for the lifting surface, including the winglet if applicable.

        Returns:
            list: A list of airfoil profiles.
        """
        profiles = [self.sections[0].root_airfoil]
        for i in range(self.num_sections):
            profiles.append(self.sections[i].tip_airfoil)
        if self.has_winglet:
            profiles.append(self.winglet.tip_airfoil)
        return profiles

    @Attribute
    def semi_span(self):
        """
        Computes the semi-span of the lifting surface.

        Returns:
            float: The semi-span of the lifting surface.
        """
        return abs(self.sections[-1].tip_airfoil.position.y - self.sections[0].root_airfoil.position.y)
    
    @Attribute
    def half_area(self):
        """
        Computes the half area of the lifting surface.

        Returns:
            float: The half area of the surface.
        """
        return sum([sec.area for sec in self.sections])
    
    @Attribute
    def wing_half_aspect_ratio(self):
        """
        Computes the half aspect ratio of the wing.

        Returns:
            float: The half aspect ratio of the wing.
        """
        return self.semi_span * self.semi_span / self.half_area
    
    @Attribute
    def mean_aerodynamic_chord(self):
        """
        Computes the mean aerodynamic chord (MAC) for the entire lifting surface.

        Returns:
            float: The mean aerodynamic chord.
        """
        return sum([sec.section_mean_aerodynamic_chord * sec.span for sec in self.sections]) / self.semi_span

    @Attribute
    def q3d_planform_geom(self) -> matlab.double:
        """
        Returns the Q3D planform geometry for the lifting surface.

        Returns:
            matlab.double: The Q3D planform geometry.
        """
        planform_geom = []
        for i in range(len(self.sections) + 1):
            x, y, z = self.profiles[i].position.location
            twist = 0.0 if i == 0 else self.sections[i-1].twist
            planform_geom.append([x, y, z, self.profiles[i].chord, twist])
        return matlab.double(planform_geom)
    
    @Attribute
    def q3d_cst_airfoils(self) -> matlab.double:
        """
        Returns the Q3D CST airfoils for the lifting surface.

        Returns:
            matlab.double: The Q3D CST airfoils.
        """
        cst_coeffs = []
        for i in range(len(self.sections) + 1):
            cst_coeffs.append(self.profiles[i].cst_coeff_u + self.profiles[i].cst_coeff_l)
        cst_coeffs = matlab.double(cst_coeffs)
        return cst_coeffs
    
    @Attribute
    def q3d_eta_airfoils(self) -> matlab.double:
        """
        Returns the Q3D eta values for the airfoils along the span of the lifting surface.

        Returns:
            matlab.double: The Q3D eta values for the airfoils.
        """
        eta = [[0.0]]
        for i in range(1, len(self.sections) + 1):
            eta.append([self.profiles[i].position.location.y / self.semi_span])
        return matlab.double(eta)
    
    @Attribute
    def sweep_4c(self):
        """
        Calculates the sweep angle at the 4c (quarter chord) of the lifting surface.

        Returns:
            float: The sweep angle at the quarter chord.
        """
        root_4c = self.profiles[0].position.x + self.profiles[0].chord / 4
        if self.has_winglet:
            tip_4c = self.profiles[-2].position.x + self.profiles[-2].chord / 4
        else:
            tip_4c = self.profiles[-1].position.x + self.profiles[-1].chord / 4
        return np.rad2deg(np.arctan((tip_4c - root_4c) / self.semi_span))
    
    @Attribute
    def x_LEMAC(self):
        """
        Calculates the location of the leading edge of the mean aerodynamic chord (MAC) 
        for the entire lifting surface.

        Returns:
            float: The location of the leading edge of the MAC.
        """
        for sec in self.sections:
            if self.mean_aerodynamic_chord <= sec.root_chord and self.mean_aerodynamic_chord >= sec.tip_chord:
                MAC_local_taper = self.mean_aerodynamic_chord / sec.root_chord
                section_eta_MAC = sec.taper_ratio / MAC_local_taper
                x_LEMAC = sec.root_airfoil.position.x + (sec.tip_airfoil.position.x - sec.root_airfoil.position.x) * section_eta_MAC
            else:
                continue
        return x_LEMAC
    
    @Attribute
    def x_ac(self):
        """
        Returns the aerodynamic center (AC) location for the lifting surface.

        Returns:
            float: The location of the aerodynamic center.
        """
        return self.x_LEMAC + 0.25 * self.mean_aerodynamic_chord

    @Part
    def frame(self):
        """
        Returns the reference frame for the lifting surface, used for positioning and transformation.

        Returns:
            Frame: The reference frame for the lifting surface.
        """
        return Frame(pos=self.position, hidden=False)