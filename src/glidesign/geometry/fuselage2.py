# Python native imports
import math


# Python third party imports


# ParaPy imports
from parapy.geom import GeomBase, Circle, BSplineCurve, Point, PointCloud, Vector, translate, rotate, MirroredShape, LoftedSurface, \
    Position, LoftedShell, FittedCurve, Revolution, RevolvedSurface
from parapy.core import Input, Attribute, Part, child
from parapy.core.validate import OneOf, Range, GreaterThan, Positive
from wx.lib.floatcanvas.FCObjects import PointSet


# Custom imports
#from .airfoil import Airfoil

#--------------------------------------
#Version under construction:
#NASA Contractor Report 3970
#Design of Fuselage Shapes
#for Natural Laminar Flow
#--------------------------------------

class GliderFuselage(GeomBase):

    name: str = Input()
    color: str = Input('white')

    # Design parameters
    xm = Input(0.35, validator = Range(0, 0.6))                     # X_m / L
    k1 = Input(0.17, validator = Positive)                                      # curvature at X_m
    rn = Input(0.349, validator = Positive)                                     # radius of curvature at nose (non-dimensional)
    ri = Input(0.469, validator = Range(0, 1))                       # profile radius at X_i
    si = Input(1.666, validator = Positive)                                     # profile slope at X_i
    xi = Input(0.855, validator = Range(0.3, 1))                     # location of inflection point (non-dimensional)
    phi = Input(6.32, validator = Range(6.3, 6.4))                    # trailing-edge angle in radians
    L = Input(7)                                                               # fuselage length in meter
    D = Input(1.0)                                                              # max diameter

    #TODO: make this fr (L/D) iterative with the analysis:
    fr = 10 #derived(lambda self: self.L / self.D)  #Lift-to-Drag ratio (from analysis)

    # --------- Forebody ---------
    def f1(self, x):
        return -2 * x * (x - 1) ** 3

    def f2(self, x):
        return -x ** 2 * (x - 1) ** 2

    def g(self, x):
        return x ** 2 * (3 * x ** 2 - 8 * x + 6)

    def r_forebody(self, x):
        x_hat = x / (self.xm * self.L) #Location along forebody (relative)
        value_fore = (1 / (2 * self.fr)) * (self.rn * self.f1(x_hat) + self.k1 * self.f2(x_hat) + self.g(x_hat))
        return math.sqrt(value_fore)

    # --------- Midbody ---------
    def f1_mid(self, x):
        return -0.5 * x ** 3 * (x - 1) ** 2

    def f2_mid(self, x):
        return x - x ** 3 * (3 * x ** 2 - 8 * x + 6)

    def g_mid(self, x):
        return x ** 3 * (6 * x ** 2 - 15 * x + 10)

    def r_midbody(self, x):
        x_hat = (self.xi * self.L - x) / ((self.xi - self.xm) * self.L) #Location along midbody (1 if at xm, 0 if at xi)
        k1m = ((self.xi / self.xm) - 1) ** 2 * self.k1 / (1 - self.ri)
        value_mid = (1 / (2 * self.fr)) * (
                self.ri + (1 - self.ri) * (k1m * self.f1_mid(x_hat) + self.si * self.f2_mid(x_hat) + self.g_mid(x_hat))
        )
        return value_mid

    # --------- Afterbody ---------
    @Attribute
    def s_ia(self):
        return ((1 - self.ri) * (1 - self.xi) * self.si) / ((self.xi - self.xm) * self.ri)

    @Attribute
    def s_il(self):
        return (2 * self.fr * (self.xi - self.xm) * math.tan(self.phi)) / (1 - self.ri)

    def r_afterbody(self, x):
        x_hat = x / self.L
        term1 = self.s_il * x_hat * (1 - x_hat ** 3)
        term2 = -self.s_ia * x_hat ** 2 * (2 * x_hat - 3) * (x_hat - 1)
        term3 = x_hat ** 2 * (3 * x_hat ** 2 - 8 * x_hat + 6)
        return term1 + term2 + term3

    # --------- Full Profile Points ---------
    @Attribute
    def profile_points(self):
        def get_r(x):
            if 0 <= x <= self.xm * self.L:
                return self.r_forebody(x)
            elif self.xm * self.L < x <= self.xi * self.L:
                return self.r_midbody(x)
            elif self.xi * self.L < x <= self.L:
                return self.r_afterbody(x)

        num_points = int(self.L * 100) #100 points per meter (1 each cm)
        dx = self.L / (num_points - 1)
        pts = [Point(x=x, y=0, z=(get_r(x) * self.D / 2)) for x in [i * dx for i in range(num_points)]]
        return pts

    @Part
    def point_line(self):
        return PointCloud(self.profile_points)

    @Part
    def fuselage_outerline(self):
        return FittedCurve(points=self.profile_points,
                           color='black'
        )

    @Part
    def fuselage_surface(self):
        return RevolvedSurface(basis_curve=self.fuselage_outerline,

                               color=self.color
        )




if __name__ == '__main__':
    from parapy.gui import display
    fus = GliderFuselage(name="fuselage")
    print(fus.profile_points)
    display(fus)
