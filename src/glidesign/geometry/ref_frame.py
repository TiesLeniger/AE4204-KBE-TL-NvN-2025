from parapy.geom import GeomBase, XOY, LineSegment, translate, Position
from parapy.core import Input, Attribute, Part, child

class Frame(GeomBase):
    """
    A class that represents a coordinate frame in 3D space, positioned at the origin 
    (default XOY plane). It has three vectors (X, Y, Z) that define the frame's axes, 
    and it visualizes these axes using line segments.

    Attributes:
        pos (Position): The position of the frame (default is the XOY plane).
        colors (list): A list of colors for visualizing the X, Y, and Z vectors.
        vectors (list): A list of vectors representing the X, Y, and Z directions of the frame.
    """
    
    # Position of the frame, defaulting to the XOY plane
    pos: Position = Input(XOY)

    @Attribute
    def colors(self):
        """
        Returns a list of colors for visualizing the X, Y, and Z vectors.

        Returns:
            list: A list of strings representing the colors ("red", "green", "blue") for the vectors.
        """
        return ["red", "green", "blue"]

    @Attribute
    def vectors(self):
        """
        Returns the vectors representing the X, Y, and Z axes of the frame.

        Returns:
            list: A list of vectors representing the X, Y, and Z directions of the frame.
        """
        return [self.pos.Vx, self.pos.Vy, self.pos.Vz]

    @Part
    def vector(self):
        """
        Generates a visual representation of the three frame vectors (X, Y, Z) as line segments.

        Returns:
            LineSegment: The line segment representing one of the three frame vectors. 
            The color and direction of the vector depend on the child index.
        """
        return LineSegment(
            quantify=3,  # Create three line segments (one for each vector)
            start=self.pos.location,  # Start point is the position of the frame
            end=translate(self.pos.location, self.vectors[child.index], 0.3),  # End point is translated in the direction of the vector
            color=self.colors[child.index],  # Color of the vector (red, green, blue)
            line_thickness=2  # Thickness of the line for the vector
        )