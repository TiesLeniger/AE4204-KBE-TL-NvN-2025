from parapy.core import Base, Input, Attribute
from parapy.core.validate import GreaterThanOrEqualTo, Range
from math import sqrt

# Define constants
P0 = 101325.0                           # [Pa], mean sea level pressure
RHO0 = 1.225                            # [kg/m^3], mean sea level density
T0 = 288.15                             # [K], mean sea level temperature
GAMMA = 1.4                             # [-], ratio of specific heats for air
GAS_CONST = 287.05                      # [J/kg K], gas constant of air
G0 = 9.80665                            # [m/s^2], gravitational acceleration
LAPSE_0 = -0.0065                       # [m^-1], lapse rate in the troposphere
MU_AIR = 1.789e-5                       # [Pa s], absolute viscosity of air

class IntStdAtmosphere(Base): 
    # Input height
    h = Input(0.0, validator = Range(0.0, 11000.0))

    @Attribute
    def temperature(self):
        return T0 + LAPSE_0 * self.h
    
    @Attribute
    def pressure(self):
        return P0*((1+LAPSE_0*self.h/T0)**(-G0/(GAS_CONST*LAPSE_0)))
    
    @Attribute
    def density(self):
        return self.pressure/(GAS_CONST * self.temperature)
    
    @Attribute
    def speed_of_sound(self):
        return sqrt(GAMMA * GAS_CONST * self.temperature)
    
    @Attribute
    def dynamic_viscosity(self):
        return MU_AIR
    
    @Attribute
    def kinematic_viscosity(self):
        return self.dynamic_viscosity / self.density
    

    
