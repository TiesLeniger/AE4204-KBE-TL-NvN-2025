from parapy.core import Base, Input, Attribute, Part
from parapy.core.validate import GreaterThan, GreaterThanOrEqualTo, LessThan
from ..core import IntStdAtmosphere

class Q3DData(Base):
    """
    A class representing the input data for a Q3D aerodynamic analysis, including
    velocity, angle of attack, altitude, and mean aerodynamic chord. It also computes
    key aerodynamic quantities such as Mach number, Reynolds number, and air density.
    
    Attributes:
        velocity (float): The velocity at which the Q3D analysis is run (in m/s).
        alpha (float): The angle of attack at which the Q3D analysis is run (in degrees).
        altitude (float): The altitude for the Q3D analysis (in meters).
        mean_aerodynamic_chord (float): The mean aerodynamic chord of the lifting surface (in meters).
        isa (IntStdAtmosphere): An instance of the IntStdAtmosphere class to provide atmospheric properties.
    """
    
    # Velocity at which the Q3D analysis is run, validated to be greater than 0
    velocity = Input(30.0, validator=GreaterThan(0.0))  # [m/s]
    
    # Angle of attack at which the Q3D analysis is run, validated to be less than 14 degrees
    alpha = Input(2.0, validator=LessThan(14.0))  # [deg]
    
    # Altitude for Q3D analysis, validated to be greater than 0 meters
    altitude = Input(500.0, validator=GreaterThan(0.0))  # [m]
    
    # Mean aerodynamic chord of the lifting surface, used for Reynolds number calculation
    mean_aerodynamic_chord = Input(0.5, validator=GreaterThan(0.0))  # [m]

    @Part
    def isa(self):
        """
        Creates an instance of the IntStdAtmosphere class with the given altitude.
        
        Returns:
            IntStdAtmosphere: An instance of the IntStdAtmosphere class that provides 
                               standard atmospheric properties based on the altitude.
        """
        return IntStdAtmosphere(
            h=self.altitude  # Set the altitude for the ISA model
        )
    
    @Attribute
    def mach_number(self):
        """
        Calculates the Mach number for the given velocity and the speed of sound 
        at the current altitude (as provided by the ISA model).
        
        The Mach number is a dimensionless quantity that represents the ratio of 
        the velocity of the object to the speed of sound in the surrounding medium.
        
        Returns:
            float: The Mach number, which helps to determine compressibility effects.
        """
        return self.velocity / self.isa.speed_of_sound
    
    @Attribute
    def reynolds_number(self):
        """
        Calculates the Reynolds number for the given velocity, air density, and mean 
        aerodynamic chord length. The Reynolds number is a dimensionless quantity used 
        in viscous flow analysis.
        
        The formula is: Re = (density * velocity * chord) / dynamic_viscosity
        
        Returns:
            float: The Reynolds number for the Q3D analysis.
        """
        return (self.isa.density * self.velocity * self.mean_aerodynamic_chord) / self.isa.dynamic_viscosity
    
    @Attribute
    def density(self):
        """
        Returns the air density at the given altitude, as provided by the ISA model.
        
        Returns:
            float: The air density (kg/m^3) at the specified altitude.
        """
        return self.isa.density  # Air density for Q3D analysis [kg/m^3]