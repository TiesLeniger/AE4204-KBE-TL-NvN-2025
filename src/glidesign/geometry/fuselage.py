# Python native imports
import math
from math import degrees

# Python third-party imports

# ParaPy imports
from parapy.geom import *
from parapy.geom.future import Common
from parapy.core import Input, Attribute, Part
from parapy.core.validate import Range, GE

# Custom imports

class GliderFuselage(GeomBase):
    """
    A class to model the fuselage of a glider. The fuselage is defined using different 
    design parameters and is divided into forebody, midbody, and afterbody regions. The 
    class calculates the corresponding profiles and generates a 3D representation of the 
    fuselage and its components.

    Attributes:
        name (str): The name of the fuselage.
        color (str): The color of the fuselage.
        xm (float): Non-dimensional location of the fuselage curvature maximum.
        k1 (float): Curvature at the location xm.
        rn (float): Non-dimensional radius of curvature at the nose.
        ri (float): Profile radius at location Xi.
        si (float): Profile slope at location Xi.
        xi (float): Non-dimensional location of the inflection point.
        t (float): Trailing gap.
        L (float): Fuselage length (in meters).
        D (float): Maximum fuselage diameter.
        mesh_deflection (float): Mesh deflection for profile accuracy.
        wing_position (tuple): Position of the wing relative to the fuselage.
        fr (float): Finesse ratio of the fuselage.
        profile_points (list): A list of points representing the fuselage profile.
        point_line (PointCloud): 3D points representing the fuselage profile line.
        fuselage_outerline (InterpolatedCurve): The outer boundary curve of the fuselage.
        fuselage_solid (RevolvedSolid): The 3D solid representation of the fuselage.
        canopy_position (Point): The position of the canopy relative to the fuselage.
        canopy_ellipse (Wire): The 2D elliptical profile of the canopy.
        canopy_spheroid (RevolvedSolid): The 3D spheroid solid for the canopy.
        canopy_intersection (IntersectedShapes): The intersection of the fuselage and canopy.
    """

    # Input design parameters
    name: str = Input()
    color: str = Input('white')

    xm = Input(0.22, validator=Range(0, 0.6))          # X_m / L
    k1 = Input(0, validator=GE(0))                      # Curvature at X_m
    rn = Input(0.7, validator=GE(0))                    # Radius of curvature at nose (non-dimensional)
    ri = Input(0.8)                                     # Profile radius at X_i
    si = Input(1.3, validator=GE(0))                    # Profile slope at X_i
    xi = Input(0.32, validator=Range(0.2, 1))           # Location of inflection point (non-dimensional)
    t = Input(0.2, validator=Range(0.05, 0.4))          # Trailing gap?
    L = Input()                                          # Fuselage length in meters
    D = Input()                                          # Max diameter
    mesh_deflection = Input()
    wing_position = Input()

    @Attribute
    def fr(self):
        """
        Calculates the finesse ratio of the fuselage based on the length (L) and diameter (D).

        Returns:
            float: The finesse ratio.
        """
        return self.L / self.D  # Finesse ratio

    # --------- Forebody ---------
    def f1(self, x):
        """
        Function for computing the forebody profile based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated value for the forebody profile.
        """
        return -2 * x * (x - 1) ** 3

    def f2(self, x):
        """
        Another function for computing the forebody profile based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated value for the forebody profile.
        """
        return -x ** 2 * (x - 1) ** 2

    def g(self, x):
        """
        A helper function for computing the curvature at the forebody based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated curvature at the forebody.
        """
        return x ** 2 * (3 * x ** 2 - 8 * x + 6)

    def r_forebody(self, x):
        """
        Calculates the forebody radius at a given position 'x' based on the fuselage design parameters.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The forebody radius.
        """
        x_hat = x / self.xm
        value_fore = (1 / (2 * self.fr)) * math.sqrt(self.rn * self.f1(x_hat) + self.k1 * self.f2(x_hat) + self.g(x_hat))
        return value_fore

    # --------- Midbody ---------
    def f1_mid(self, x):
        """
        Function for computing the midbody profile based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated value for the midbody profile.
        """
        return -0.5 * x ** 3 * (x - 1) ** 2

    def f2_mid(self, x):
        """
        Another function for computing the midbody profile based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated value for the midbody profile.
        """
        return x - x ** 3 * (3 * x ** 2 - 8 * x + 6)

    def g_mid(self, x):
        """
        A helper function for computing the curvature at the midbody based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated curvature at the midbody.
        """
        return x ** 3 * (6 * x ** 2 - 15 * x + 10)

    def r_midbody(self, x):
        """
        Calculates the midbody radius at a given position 'x' based on the fuselage design parameters.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The midbody radius.
        """
        x_hat = (self.xi - x) / (self.xi - self.xm)
        k1m = ((self.xi / self.xm) - 1) ** 2 * self.k1 / (1 - self.ri)
        value_mid = (1 / (2 * self.fr)) * (self.ri + (1 - self.ri) * (k1m * self.f1_mid(x_hat) + self.si * self.f2_mid(x_hat) + self.g_mid(x_hat)))
        return value_mid

    # --------- Afterbody ---------
    @Attribute
    def s_ia(self):
        """
        Computes a parameter used for the afterbody profile based on the fuselage design parameters.

        Returns:
            float: The computed s_ia value for the afterbody.
        """
        return ((1 - self.ri) * (1 - self.xi) * self.si) / ((self.xi - self.xm) * self.ri)

    def f1_aft(self, x):
        """
        Function for computing the afterbody profile based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated value for the afterbody profile.
        """
        return 1 - x ** 3 * (6 * x ** 2 - 15 * x + 10)

    def f2_aft(self, x):
        """
        Another function for computing the afterbody profile based on the parameter 'x'.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The calculated value for the afterbody profile.
        """
        return -x ** 3 * (3 * x ** 2 - 7 * x + 4)

    def r_afterbody(self, x):
        """
        Calculates the afterbody radius at a given position 'x' based on the fuselage design parameters.

        Args:
            x (float): A non-dimensional parameter representing the position along the fuselage.

        Returns:
            float: The afterbody radius.
        """
        x_hat = (1 - x) / (1 - self.xi)
        value_aft = (self.ri / (2 * self.fr)) * (1 + ((self.t / self.ri) - 1) * self.f1_aft(x_hat) + self.s_ia * self.f2_aft(x_hat))
        return value_aft

    # --------- Full Profile Points ---------
    @Attribute
    def profile_points(self):
        """
        Computes the full set of profile points for the fuselage, combining the forebody, midbody, 
        and afterbody radii based on the position 'x_hat'.

        Returns:
            list: A list of points representing the fuselage profile.
        """
        def get_r(x_hat):
            if 0 <= x_hat <= self.xm:
                return self.r_forebody(x_hat)
            elif self.xm < x_hat <= self.xi:
                return self.r_midbody(x_hat)
            elif self.xi < x_hat <= 1:
                return self.r_afterbody(x_hat)

        num_points = int(self.L * 10)  # 10 points per meter (1 each 10 cm)
        dx_hat = 1.0 / (num_points - 1)

        pts = [Point(x=x_hat * self.L,
                     y=0,
                     z=get_r(x_hat) * self.L)
               for x_hat in [i * dx_hat for i in range(num_points)]]
        return pts

    @Attribute
    def point_line(self):
        """
        Returns a point cloud representing the fuselage profile.

        Returns:
            PointCloud: A point cloud of the fuselage profile.
        """
        return PointCloud(self.profile_points)

    @Attribute
    def fuselage_outerline(self):
        """
        Returns the outer boundary curve of the fuselage as an interpolated curve.

        Returns:
            InterpolatedCurve: The outer boundary curve of the fuselage.
        """
        return InterpolatedCurve(points=self.profile_points,
                                 color='black')

    @Part
    def fuselage_solid(self):
        """
        Returns the 3D solid representation of the fuselage using the revolved outerline.

        Returns:
            RevolvedSolid: The 3D solid representation of the fuselage.
        """
        return RevolvedSolid(built_from=self.fuselage_outerline,
                              center=self.position,
                              direction=Vector(1, 0, 0),
                              color=self.color)

    # Canopy
    @Attribute
    def canopy_position(self):
        """
        Returns the position of the canopy relative to the fuselage.

        Returns:
            Point: The position of the canopy.
        """
        return translate(self.position,
                         'x', (self.xm - 0.08) * self.L,
                         'z', self.D)

    @Attribute
    def canopy_ellipse(self):
        """
        Returns the 2D elliptical profile of the canopy.

        Returns:
            Wire: The 2D elliptical profile of the canopy.
        """
        return Wire([Ellipse(position=self.canopy_position,
                             major_radius=self.wing_position[0] / 2,  # lengthwise (X-axis)
                             minor_radius=self.D)], hidden=True)  # heightwise (Z-axis)

    @Attribute
    def canopy_spheroid(self):
        """
        Returns the 3D spheroid representation of the canopy.

        Returns:
            RevolvedSolid: The 3D spheroid solid for the canopy.
        """
        return RevolvedSolid(center=self.canopy_position,
                             built_from=self.canopy_ellipse,
                             direction=Vector(1, 0, 0),
                             color="RED",
                             hidden=True)

    @Part
    def canopy_intersection(self):
        """
        Returns the intersection of the fuselage and the canopy.

        Returns:
            IntersectedShapes: The intersection of the fuselage and canopy shapes.
        """
        return IntersectedShapes(self.fuselage_solid, self.canopy_spheroid, color='red')