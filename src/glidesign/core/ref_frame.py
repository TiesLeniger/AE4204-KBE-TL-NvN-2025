from parapy.geom import GeomBase, XOY, LineSegment, translate, Position
from parapy.core import Input, Attribute, Part, child

class Frame(GeomBase):
    pos: Position = Input(XOY)

    @Attribute
    def colors(self):
        return ["red", "green", "blue"]

    @Attribute
    def vectors(self):
        return [self.pos.Vx, self.pos.Vy, self.pos.Vz]

    @Part
    def vector(self):
        return LineSegment(quantify=3,
                           start=self.pos.location,
                           end=(translate(self.pos.location, self.vectors[child.index], 0.3)),
                           color=self.colors[child.index],
                           line_thickness=2
                           )