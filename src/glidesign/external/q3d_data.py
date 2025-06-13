from parapy.core import Base, Input, Attribute, Part
from parapy.core.validate import GreaterThan, GreaterThanOrEqualTo, LessThan
from ..core import IntStdAtmosphere

class Q3DData(Base):

    velocity = Input(30.0, validator = GreaterThan(0.0))                # Velocity at which Q3D analysis is run [m/s]
    alpha = Input(2.0, validator = LessThan(14.0))                      # Angle of attack at which Q3D is run [deg]
    altitude = Input(500.0, validator = GreaterThan(0.0))               # Altitude for Q3D [m]
    mean_aerodynamic_chord = Input(0.5, validator = GreaterThan(0.0))   # Mean aerodynamic chord of lifting surface, used for Re

    @Part
    def isa(self):
        return IntStdAtmosphere(
            h = self.altitude
        )
    
    @Attribute
    def mach_number(self):                                          # Mach number for adding compressibility effects to Q3D [-]
        return self.velocity / self.isa.speed_of_sound
    
    @Attribute
    def reynolds_number(self):                                      # Reynolds number for viscous analysis of Q3D [-]
        return (self.isa.density * self.velocity * self.mean_aerodynamic_chord) / self.isa.dynamic_viscosity
    
    @Attribute
    def density(self):
        return self.isa.density                                     # Air density for Q3D analysis [kg/m^3]