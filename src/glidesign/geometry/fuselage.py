# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, Circle, BSplineCurve, Point, Vector, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part, child
from parapy.core.validate import OneOf, Range, GreaterThan

# Custom imports
#from .airfoil import Airfoil

class GliderFuselage(GeomBase):

    name: str = Input()
    color: str = Input('white')

    length: float = Input()
    max_diameter: float = Input()

    relative_section_diameter: list[float] = Input([10, 50, 80, 95, 100, 75, 45, 25, 15, 10], validator = Range(1, 100)) #Percentage of maximum fuselage diameter

    @Attribute
    def section_diameter(self) -> list[float]:
        return [i * self.max_diameter/100. for i in self.relative_section_diameter]

    @Attribute
    def section_length(self) -> float:
        return self.length / (len(self.relative_section_diameter) - 1)

    @Attribute
    def spine_points(self):
        return [
            Point(0.2432, 0, 0.0021),
            Point(0.5012, 0, 0.0144),
            Point(0.7532, 0, 0.0144),
            Point(1.0053, 0, 0.0433),
            Point(1.2485, 0, 0.0680),
            Point(1.4739, 0, 0.0948),
            Point(1.6904, 0, 0.1196),
            Point(1.8831, 0, 0.1402),
            Point(2.1263, 0, 0.1629),
            Point(2.3280, 0, 0.1649),
            Point(2.5177, 0, 0.1629),
            Point(2.7402, 0, 0.1608),
            Point(2.9033, 0, 0.1526),
            Point(3.1109, 0, 0.1402),
            Point(3.3184, 0, 0.1196),
            Point(3.5468, 0, 0.1052),
            Point(3.7870, 0, 0.1010),
            Point(4.0391, 0, 0.0887),
            Point(4.2823, 0, 0.0845),
            Point(4.5432, 0, 0.0804),
            Point(4.8072, 0, 0.0722),
            Point(5.1037, 0, 0.0680),
            Point(5.3884, 0, 0.0701),
            Point(5.7205, 0, 0.0639),
            Point(6.0112, 0, 0.0619),
            Point(6.2959, 0, 0.0639) #Points retreived from approximately the AS33 fuselage
        ]

    @Attribute
    def fuselage_spine(self):
        return BSplineCurve(points=self.spine_points,
                                 color='black'
                                )

    @Attribute
    def x_equidistant_points_bspline(self):
        num_samples = 1000
        t_samples = np.linspace(0, 1, num_samples)
        points = [self.bspline_curve.position_at(t).point for t in t_samples]
        x_vals = [pt.x for pt in points]

        x_start = 0
        x_end = self.length
        x_targets = np.linspace(x_start, x_end, len(self.relative_section_diameter))

        t_targets = np.interp(x_targets, x_vals, t_samples)

        return [self.bspline_curve.position_at(t) for t in t_targets]

    @Part
    def fuselage_profiles(self):
        return Circle(quantify=len(self.sections),
                      color = 'white',
                      radius = self.section_diameter[child.index]/2, #Convert diameter into radius
                      position = translate(self.position.rotate90('y'),
                                           self.x_equidistant_points_bspline[child.index])
                      )


