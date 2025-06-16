# Python native imports
import os

# Python third-party package imports
import numpy as np

# Parapy imports
from parapy.geom import FittedCurve, Point
from parapy.core import Attribute, Part, Input
from parapy.core.validate import OneOf, Range, GreaterThan
from parapy.gui import display
from kbeutils.data import airfoils
from kbeutils.geom import Naca4AirfoilCurve, Naca5AirfoilCurve
import kbeutils.avl as avl

# Self-built imports
from ..core import airfoil_found, NACA4_PATTERN, NACA5_PATTERN
from ..core.cst_curves import fit_cst_airfoil
from .ref_frame import Frame

class Airfoil(FittedCurve):
    """
    A class representing an airfoil, which can either be defined by a NACA 4-digit, 
    NACA 5-digit designation or a custom airfoil file. The class computes the airfoil 
    coordinates, CST coefficients, and provides methods for interfacing with aerodynamic 
    analysis tools.

    Attributes:
        airfoil_name (str): The name of the airfoil (NACA 4, NACA 5, or custom).
        chord (float): The chord length of the airfoil (in meters).
        cst_poly_order (int): The polynomial order for CST fitting.
        num_points (int): The number of points to generate for the airfoil coordinates.
        mesh_deflection (float): The mesh deflection used for airfoil geometry generation.
        coords (list): The coordinates of the airfoil (x, z).
        points (list[Point]): The points of the airfoil in 3D space.
        cst_coefficients (tuple): The CST coefficients for the upper and lower surfaces of the airfoil.
        cst_coeff_u (list): The CST coefficients for the upper surface.
        cst_coeff_l (list): The CST coefficients for the lower surface.
        z_te_u (float): The trailing edge location for the upper surface.
        z_te_l (float): The trailing edge location for the lower surface.
        frame (Frame): A reference frame for the airfoil position.
        avl_section (avl.SectionFromCurve): The AVL section representation of the airfoil.
    """

    # Input values for the airfoil
    airfoil_name: str = Input("nlf1-0015", validator=airfoil_found)
    chord: float = Input(1.0, validator=GreaterThan(0.0))
    cst_poly_order: int = Input(5, validator=GreaterThan(0))
    num_points: int = Input(200, validator=GreaterThan(0))
    mesh_deflection: float = Input(1e-5, validator=GreaterThan(0.0))

    @Attribute
    def coords(self):
        """
        Retrieves the coordinates of the airfoil either from the NACA 4-digit or NACA 5-digit 
        airfoil definition or from a custom airfoil file. The coordinates are returned as normalized 
        x and z values relative to the chord length.

        Returns:
            list: A list containing the x and z coordinates of the airfoil.
        """
        # Check if airfoil is NACA 4-digit or 5-digit, otherwise, fetch custom coordinates
        if NACA4_PATTERN.match(self.airfoil_name):
            digits = NACA4_PATTERN.match(self.airfoil_name).group(1)
            naca4curve = Naca4AirfoilCurve(designation=digits, mesh_deflection=1e-6, hidden=True)
            coordinates = np.array(np.array(naca4curve.coordinates))
            assert coordinates.shape == (naca4curve.n_points, 3), f"Expected shape ({naca4curve.n_points}, 3), got {coordinates.shape}"
            coordinates = np.round(coordinates / naca4curve.chord_length, decimals=5)
            return [coordinates[:, 0].tolist(), coordinates[:, 2].tolist()]
        
        elif NACA5_PATTERN.match(self.airfoil_name):
            digits = NACA5_PATTERN.match(self.airfoil_name).group(1)
            naca5curve = Naca5AirfoilCurve(designation=digits, mesh_deflection=1e-6, hidden=True)
            coordinates = np.array(np.array(naca5curve.coordinates))
            assert coordinates.shape == (naca5curve.n_points, 3), f"Expected shape ({naca5curve.n_points}, 3), got {coordinates.shape}"
            coordinates = np.round(coordinates / naca5curve.chord_length, decimals=5)
            return [coordinates[:, 0].tolist(), coordinates[:, 2].tolist()]
        
        else:
            # Load airfoil coordinates from file if custom airfoil is provided
            in_input = os.path.exists(os.path.join(os.getcwd(), "input", "airfoils", self.airfoil_name + ".dat"))
            path_to_af_file = os.path.join(os.getcwd(), "input", "airfoils", self.airfoil_name + ".dat") if in_input else os.path.join(
                airfoils.__path__[0], self.airfoil_name + ".dat")
            with open(path_to_af_file, 'r') as file:
                x_coords_lst = []
                z_coords_lst = []
                skipped_lines = 0
                for line in file:
                    try:
                        x, z = line.split(maxsplit=1)
                        x_coords_lst.append(float(x.strip()))
                        z_coords_lst.append(float(z.strip()))
                    except ValueError:
                        if skipped_lines == 0:
                            skipped_lines = 1
                            continue
                        else:
                            break
            return [x_coords_lst, z_coords_lst]

    @Attribute
    def points(self) -> list[Point]:
        """
        Returns the 3D points of the airfoil in space based on its coordinates.
        
        Returns:
            list[Point]: A list of Parapy `Point` objects representing the airfoil geometry.
        """
        return [self.position.translate("x", x * self.chord, "z", z * self.chord).location
                for x, z in zip(self.coords[0], self.coords[1])]

    @Attribute
    def cst_coefficients(self):
        """
        Returns the CST coefficients for the airfoil, based on the specified polynomial order.

        Returns:
            tuple: The CST coefficients for the upper and lower surfaces and trailing edge positions.
        """
        return fit_cst_airfoil(np.array(self.coords[0]), np.array(self.coords[1]), self.cst_poly_order)

    @Attribute
    def cst_coeff_u(self):
        """
        Returns the CST coefficients for the upper surface of the airfoil.
        
        Returns:
            list: The CST coefficients for the upper surface.
        """
        return self.cst_coefficients[0]
    
    @Attribute
    def cst_coeff_l(self):
        """
        Returns the CST coefficients for the lower surface of the airfoil.
        
        Returns:
            list: The CST coefficients for the lower surface.
        """
        return self.cst_coefficients[1]
    
    @Attribute
    def z_te_u(self):
        """
        Returns the trailing edge location for the upper surface.
        
        Returns:
            float: The trailing edge position for the upper surface.
        """
        return self.cst_coefficients[2]
    
    @Attribute
    def z_te_l(self):
        """
        Returns the trailing edge location for the lower surface.
        
        Returns:
            float: The trailing edge position for the lower surface.
        """
        return self.cst_coefficients[3]
    
    @Part
    def frame(self):
        """
        Returns a reference frame for the airfoil, used for positioning and transformation.
        
        Returns:
            Frame: The reference frame for the airfoil.
        """
        return Frame(pos=self.position, hidden=False)
    
    @Part
    def avl_section(self):
        """
        Returns the AVL section representation of the airfoil, used for aerodynamic analysis.
        
        Returns:
            avl.SectionFromCurve: The AVL section based on the airfoil geometry.
        """
        return avl.SectionFromCurve(
            curve_in=self  # Pass the current airfoil as the curve input for AVL analysis
        )