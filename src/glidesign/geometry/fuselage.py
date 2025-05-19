# Python native imports

# Python third party imports
import numpy as np

# ParaPy imports
from parapy.geom import GeomBase, Circle, InterpolatedCurve, Point, Vector, translate, rotate, MirroredShape
from parapy.core import Input, Attribute, Part, child
from parapy.core.validate import OneOf, Range, GreaterThan

# Custom imports
from ..core import airfoil_found
from .airfoil import Airfoil

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
            Point(6.2959, 0, 0.0639)
        ]

    @Attribute
    def fuselage_spine(self):
        return InterpolatedCurve(points=self.spine_points,
                                 initial_tangent=Vector(1, 0, 0),
                                 final_tangent=Vector(1, 0, 0),
                                 color='black'
                                )

    @Part
    def fuselage_profiles(self):
        return Circle(quantify=len(self.sections),
                      color = 'white',
                      radius = self.section_diameter[child.index]/2 #Convert diameter into radius
                      position =
                      )


