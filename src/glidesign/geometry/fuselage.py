# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, Circle, BSplineCurve, Point, Vector, translate, rotate, MirroredShape, LoftedSurface, \
    Position, LoftedShell
from parapy.core import Input, Attribute, Part, child
from parapy.core.validate import OneOf, Range, GreaterThan

# Custom imports
#from .airfoil import Airfoil

class GliderFuselage(GeomBase):

    name: str = Input()
    color: str = Input('white')

    length: float = Input()
    max_diameter: float = Input()

    relative_section_diameter: list[float] = Input([5, 40, 60, 80, 90, 100, 100, 100, 100, 90, 80, 70, 60, 50, 40, 35, 33, 29, 28, 27, 26, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25]) #Percentage of maximum fuselage diameter

    @Attribute
    def section_diameter(self) -> list[float]:
        return [i * self.max_diameter/100. for i in self.relative_section_diameter]

    @Attribute
    def section_length(self) -> float:
        return self.length / (len(self.relative_section_diameter) - 1)

    @Attribute
    def spine_points(self):
        return [
            Point(0.0, 0, 0.0),               # nose
            Point(1.5, 0, 0.2),               # cockpit rise
            Point(2.3, 0, 0.4),
            Point(3.0, 0, 0.6),               # max diameter
            Point(4.5, 0, 0.6),               # beginning taper
            Point(6.0, 0, 0.6),               # tail boom
            Point(self.length, 0, 0.6)                # tail tip
        ]

    @Part
    def fuselage_spine(self):
        return BSplineCurve(self.spine_points,
                            color='black'
                            )

    @Attribute
    def x_equidistant_points_bspline(self):
        num_samples = 1000
        t_samples = np.linspace(0, 1, num_samples)
        points = [self.fuselage_spine.point_at_parameter(t) for t in t_samples]
        x_vals = [pt.x for pt in points]

        x_start = 0
        x_end = self.length
        x_targets = np.linspace(x_start, x_end, len(self.relative_section_diameter))

        t_targets = np.interp(x_targets, x_vals, t_samples)

        return [self.fuselage_spine.point_at_parameter(t) for t in t_targets]

    @Part
    def fuselage_profiles(self):
        return Circle(quantify=len(self.relative_section_diameter),
                      color = 'black',
                      radius = self.section_diameter[child.index]/2, #Convert diameter into radius
                      position=translate(self.position,
                                         'x', self.x_equidistant_points_bspline[child.index].x,
                                         'y', self.x_equidistant_points_bspline[child.index].y,
                                         'z', self.x_equidistant_points_bspline[child.index].z
                                         ).rotate('y', 90, deg=True)
                      # Apply per-profile rotation
                      )

    @Part
    def fuselage_lofted(self):
        return LoftedShell(profiles = self.fuselage_profiles,
                           ruled=True,
                           color = self.color,
                           hidden = not(__name__ == '__main__')
        )

if __name__ == '__main__':
    from parapy.gui import display
    fus = GliderFuselage(name="fuselage", max_diameter= 0.75, length= 6.5)
    print(fus.x_equidistant_points_bspline)
    display(fus)
