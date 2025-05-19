from parapy.core import *
from parapy.geom import *
import math

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
    xm = Input(0.4)     # X_m / L
    k1 = Input(0.2)     # curvature at X_m
    rn = Input(0.05)    # radius of curvature at nose (non-dimensional)
    ri = Input(0.6)     # profile radius at X_i
    si = Input(0.5)     # profile slope at X_i
    xi = Input(0.8)     # location of inflection point (non-dimensional)
    phi = Input(math.radians(10))  # trailing-edge angle in radians
    L = Input(10.0)     # fuselage length
    D = Input(1.0)      # max diameter

    fr = 50 #derived(lambda self: self.L / self.D)  # f_r

    # --------- Forebody ---------
    @Attribute
    def F1(self):
        return -2 * (x - 1) ** 3

    @Attribute
    def F2(self):
        return -x**2 * (x - 1) ** 2

    @Attribute
    def G(self):
        return x**2 * (3 * x**2 - 8 * x + 6)


    def r_forebody(self):
        x_hat = x / (self.xm * self.L)
        val = (1 / (2 * self.fr)) * (self.rn * self.F1(x_hat) + self.k1 * self.F2(x_hat) + self.G(x_hat))
        return math.sqrt(val)

    # --------- Midbody ---------
    def F1_mid(self, x): return -0.5 * x**3 * (x - 1)**2
    def F2_mid(self, x): return -x**3 * (3 * x**2 - 8 * x + 6)
    def G_mid(self, x): return x**3 * (6 * x**2 - 15 * x + 10)

    def r_midbody(self, x):
        x_hat = (self.xi * self.L - x) / ((self.xi - self.xm) * self.L)
        k1m = ((self.xi / self.xm) - 1) ** 2 * self.k1 / (1 - self.ri)
        val = (1 / (2 * self.fr)) * (
            self.ri + (1 - self.ri) * (k1m * self.F1_mid(x_hat) + self.si * self.F2_mid(x_hat) + self.G_mid(x_hat))
        )
        return val

    # --------- Afterbody ---------
    @Attribute
    def s_ia(self):
        return ((1 - self.ri) * (1 - self.xi) * self.si) / ((self.xi - self.xm) * self.ri)

    @Attribute
    def s_iL(self):
        return 2 * self.fr * (self.xi - self.xm) * math.tan(self.phi) / (1 - self.ri)

    def r_afterbody(self, x):
        x_hat = x / self.L
        term1 = self.s_iL * x_hat * (1 - x_hat ** 3)
        term2 = -self.s_ia * x_hat ** 2 * (2 * x_hat - 3) * (x_hat - 1)
        term3 = x_hat ** 2 * (3 * x_hat ** 2 - 8 * x_hat + 6)
        return term1 + term2 + term3

    # --------- Full Profile Points ---------
    @Part
    def profile_points(self):
        def get_r(x):
            if 0 <= x <= self.xm * self.L:
                return self.r_forebody(x)
            elif self.xm * self.L < x <= self.xi * self.L:
                return self.r_midbody(x)
            else:
                return self.r_afterbody(x)

        num_points = 100
        dx = self.L / (num_points - 1)
        pts = [Point(x=x, y=0, z=get_r(x) * self.D / 2) for x in [i * dx for i in range(num_points)]]
        return pts

    @Part
    def fuselage_surface(self):
        return RevolvedSurface(
            profile=Polyline(points=self.profile_points),
            axis=Line(Point(0, 0, 0), Point(self.L, 0, 0)),
            angle=360
        )




if __name__ == '__main__':
    from parapy.gui import display
    fus = GliderFuselage(name="fuselage")
    display(fus)
