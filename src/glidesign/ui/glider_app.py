from parapy.core import *
from glidesign.geometry.glider import Glider
from glidesign.core.excel_loader import load_glider_params

class GliderApp(Base):
    file_to_load = Input("glider.xlsx")

    @Attribute
    def parameters(self):
        return load_glider_params(self.file_to_load)

    @Part
    def glider(self): 
        return Glider(
            occupants = int(self.parameters["Occupants"]),
            fai_class = str(self.parameters["FAI Class"]),
            open_class_wingspan = float(self.parameters["Open class wingspan"]),
            current_pilot_mass = float(self.parameters["Pilot mass"]),
            glider_structural_material = str(self.parameters["Structural material"])
            )