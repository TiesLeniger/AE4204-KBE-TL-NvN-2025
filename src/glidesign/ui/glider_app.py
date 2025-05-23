from parapy.core import *
from glidesign.geometry.glider import Glider
from glidesign.core.excel_loader import load_glider_params

class GliderApp(Base):
    parameters: dict = Input(load_glider_params("input/glider.xlsx"))

    @Part
    def glider(self):
        return Glider(parameters=self.parameters)