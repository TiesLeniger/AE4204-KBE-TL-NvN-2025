from parapy.core import Base, Input, Attribute
from parapy.core.validate import GreaterThanOrEqualTo, Range
from math import sqrt

# Define constants for standard atmospheric conditions
P0 = 101325.0                           # [Pa], mean sea level pressure
RHO0 = 1.225                            # [kg/m^3], mean sea level density
T0 = 288.15                             # [K], mean sea level temperature
GAMMA = 1.4                             # [-], ratio of specific heats for air
GAS_CONST = 287.05                      # [J/kg K], gas constant of air
G0 = 9.80665                            # [m/s^2], gravitational acceleration
LAPSE_0 = -0.0065                       # [m^-1], lapse rate in the troposphere
MU_AIR = 1.789e-5                       # [Pa s], absolute viscosity of air

class IntStdAtmosphere(Base):
    """
    Class that models the International Standard Atmosphere (ISA) conditions
    based on input height in the troposphere (up to 11,000 meters). 
    
    Attributes:
        h (float): The height in meters above sea level (must be between 0 and 11,000 meters).
    
    Provides the following standard atmospheric properties at the given height:
        - Temperature
        - Pressure
        - Density
        - Speed of sound
        - Dynamic viscosity
        - Kinematic viscosity
    """
    
    # Input for height (h) with a valid range from 0.0 to 11000.0 meters
    h = Input(0.0, validator=Range(0.0, 11000.0))

    @Attribute
    def temperature(self):
        """
        Calculate the temperature at the given height using the lapse rate.

        Returns:
            float: Temperature in Kelvin at the input height.
        """
        return T0 + LAPSE_0 * self.h
    
    @Attribute
    def pressure(self):
        """
        Calculate the pressure at the given height using the barometric formula.

        Returns:
            float: Pressure in Pascals at the input height.
        """
        return P0 * ((1 + LAPSE_0 * self.h / T0) ** (-G0 / (GAS_CONST * LAPSE_0)))
    
    @Attribute
    def density(self):
        """
        Calculate the air density at the given height using the ideal gas law.

        Returns:
            float: Density in kg/m^3 at the input height.
        """
        return self.pressure / (GAS_CONST * self.temperature)
    
    @Attribute
    def speed_of_sound(self):
        """
        Calculate the speed of sound at the given height using the ideal gas law.

        Returns:
            float: Speed of sound in m/s at the input height.
        """
        return sqrt(GAMMA * GAS_CONST * self.temperature)
    
    @Attribute
    def dynamic_viscosity(self):
        """
        Return the dynamic viscosity of air (constant value).

        Returns:
            float: Dynamic viscosity of air in Pa.s.
        """
        return MU_AIR
    
    @Attribute
    def kinematic_viscosity(self):
        """
        Calculate the kinematic viscosity at the given height.

        Returns:
            float: Kinematic viscosity in m^2/s at the input height.
        """
        return self.dynamic_viscosity / self.density