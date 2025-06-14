# Python native imports
import math
from math import degrees

# Python third party imports


# ParaPy imports
from parapy.geom import *
from parapy.geom.future import Common
from parapy.core import Input, Attribute, Part
from parapy.core.validate import Range, GE

# Custom imports

class GliderFuselage(GeomBase):

    name: str = Input()
    color: str = Input('white')

    # Design parameters - see FIGURE 1 - https://cafe.foundation/v2/pdf_tech/Drag.Reduction/5.AIAA-48131-445.pdf

    xm = Input(0.22, validator = Range(0, 0.6))                         # X_m / L
    k1 = Input(0, validator = GE(0))                                                # curvature at X_m
    rn = Input(0.7, validator = GE(0))                                              # radius of curvature at nose (non-dimensional)
    ri = Input(0.8)                                                                 # profile radius at X_i
    si = Input(1.3, validator = GE(0))                                              # profile slope at X_i
    xi = Input(0.32, validator = Range(0.2, 1))                         # location of inflection point (non-dimensional)
    t = Input(0.2, validator = Range(0.05, 0.4))                                                                  # Trailing gap?
    L = Input()                                 # fuselage length in meter
    D = Input()                                 # max diameter
    mesh_deflection = Input()
    wing_position = Input()

    @Attribute
    def fr(self):
        return self.L / self.D                                                      #Finesse ratio

    # --------- Forebody ---------
    def f1(self, x):
        return -2 * x * (x - 1) ** 3

    def f2(self, x):
        return -x ** 2 * (x - 1) ** 2

    def g(self, x):
        return x ** 2 * (3 * x ** 2 - 8 * x + 6)

    def r_forebody(self, x):
        x_hat = x / self.xm
        value_fore = (1 / (2 * self.fr)) * math.sqrt(self.rn * self.f1(x_hat) + self.k1 * self.f2(x_hat) + self.g(x_hat))
        return value_fore

    # --------- Midbody ---------
    def f1_mid(self, x):
        return -0.5 * x ** 3 * (x - 1) ** 2

    def f2_mid(self, x):
        return x - x ** 3 * (3 * x ** 2 - 8 * x + 6)

    def g_mid(self, x):
        return x ** 3 * (6 * x ** 2 - 15 * x + 10)

    def r_midbody(self, x):
        x_hat = (self.xi - x) / (self.xi - self.xm)
        k1m = ((self.xi / self.xm) - 1) ** 2 * self.k1 / (1 - self.ri)
        value_mid = (1 / (2 * self.fr)) * (self.ri + (1 - self.ri) * (k1m * self.f1_mid(x_hat) + self.si * self.f2_mid(x_hat) + self.g_mid(x_hat)))
        return value_mid

    # --------- Afterbody ---------
    @Attribute
    def s_ia(self):
        return ((1 - self.ri) * (1 - self.xi) * self.si) / ((self.xi - self.xm) * self.ri)

    def f1_aft(self, x):
        return 1 - x ** 3 * (6 * x ** 2 - 15 * x + 10)

    def f2_aft(self, x):
        return -x ** 3 * (3 * x ** 2 - 7 * x + 4)

    def r_afterbody(self, x):
        x_hat = (1 - x) / (1 - self.xi)
        value_aft = (self.ri / (2 * self.fr)) * (1 + ((self.t / self.ri) - 1) * self.f1_aft(x_hat) + self.s_ia * self.f2_aft(x_hat))
        return value_aft

    # --------- Full Profile Points ---------
    @Attribute
    def profile_points(self):
        def get_r(x_hat):
            if 0 <= x_hat <= self.xm:
                return self.r_forebody(x_hat)
            elif self.xm < x_hat <= self.xi:
                return self.r_midbody(x_hat)
            elif self.xi < x_hat <= 1:
                return self.r_afterbody(x_hat)

        num_points = int(self.L * 10) #10 points per meter (1 each 10 cm)
        dx_hat = 1.0 / (num_points - 1)

        pts = [Point(x=x_hat * self.L,
                     y=0,
                     z=get_r(x_hat) * self.L)
               for x_hat in [i * dx_hat for i in range(num_points)]]
        return pts

    @Attribute
    def point_line(self):
        return PointCloud(self.profile_points)

    @Attribute
    def fuselage_outerline(self):
        return InterpolatedCurve(points=self.profile_points,
                                 color='black'
        )

    @Part
    def fuselage_solid(self):
        return RevolvedSolid(built_from=self.fuselage_outerline,
                      center=self.position,
                      direction=Vector(1, 0, 0),
                      color=self.color
                      )

    #Make a canopy
    @Attribute
    def canopy_position(self):
        return translate(self.position,
                         'x', (self.xm - 0.08) * self.L,
                         'z', self.D)

    @Part
    def canopy_ellipse(self):
        return Wire([Ellipse(position = self.canopy_position,
                            major_radius=self.wing_position[0] / 2,  # lengthwise (X-axis)
                           minor_radius=self.D)],
                    hidden = True)  # heightwise (Z-axis)

    @Part
    def canopy_spheroid(self):
        return RevolvedSolid(center = self.canopy_position,
                               built_from =self.canopy_ellipse,
                               direction=Vector(1, 0, 0),
                             color = "RED",
                             hidden = True)

    @Part
    def canopy_intersection(self):
        return IntersectedShapes(self.fuselage_solid, self.canopy_spheroid, color='red')
