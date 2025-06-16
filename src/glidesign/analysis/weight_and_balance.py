# Python native imports

# Python third party imports

# ParaPy imports
from parapy.core import Base, Input, Attribute
from parapy.core.validate import Range
from parapy.core.widgets import Dropdown

# Custom imports

class WeightAndBalance(Base):
    """
    A class that calculates and manages the weight and balance of a glider, 
    ensuring that the center of gravity (CG) is within acceptable limits based 
    on the pilot's mass, structure material, and other glider components.
    
    Attributes:
        min_pilot_mass (float): The minimum allowed mass for the pilot.
        max_pilot_mass (float): The maximum allowed mass for the pilot.
        current_pilot_mass (float): The current mass of the pilot.
        glider_structure_material (str): The material used for the glider structure (e.g., "Carbon fibre", "Glass fibre").
        right_wing_cog (list): The center of gravity (CG) position of the right wing (x, y, z).
        left_wing_cog (list): The CG position of the left wing (x, y, z).
        vertical_tail_cog (list): The CG position of the vertical tail (x, y, z).
        right_hor_tail_cog (list): The CG position of the right horizontal tail (x, y, z).
        left_hor_tail_cog (list): The CG position of the left horizontal tail (x, y, z).
        fuselage_cog (list): The CG position of the fuselage (x, y, z).
        right_wing_volume (float): The volume of the right wing.
        left_wing_volume (float): The volume of the left wing.
        vertical_tail_volume (float): The volume of the vertical tail.
        right_hor_tail_volume (float): The volume of the right horizontal tail.
        left_hor_tail_volume (float): The volume of the left horizontal tail.
        fuselage_volume (float): The volume of the fuselage.
        pilot_cog (list): The center of gravity position of the pilot (x, y, z).
    """
    
    # Input values for mass and structure material
    min_pilot_mass = Input()                 # Minimum allowed mass for the pilot
    max_pilot_mass = Input()                 # Maximum allowed pilot mass
    current_pilot_mass = Input()             # Current pilot mass

    glider_structure_material = Input()      # Structure material type (e.g., "Carbon fibre", "Glass fibre")

    # Input values for center of gravity (CG) positions and volumes of various components
    right_wing_cog = Input()
    left_wing_cog = Input()
    vertical_tail_cog = Input()
    right_hor_tail_cog = Input()
    left_hor_tail_cog = Input()
    fuselage_cog = Input()

    # Input values for volumes of components
    right_wing_volume = Input()
    left_wing_volume = Input()
    vertical_tail_volume = Input()
    right_hor_tail_volume = Input()
    left_hor_tail_volume = Input()
    fuselage_volume = Input()

    pilot_cog = Input()                     # Pilot CG based on the canopy center (ellipse)

    @Attribute
    def material_density(self):
        """
        Returns the density of the material used for the glider structure.
        
        Returns:
            float: Material density (kg/m^3).
        """
        if self.glider_structure_material == "Carbon fibre":
            return 1700  # kg/m^3
        elif self.glider_structure_material == "Glass fibre":
            return 2000  # kg/m^3

    @Attribute
    def list_of_volumes(self):
        """
        Returns a list of the effective volumes for different glider components 
        based on a volumetric structure fraction.
        
        Returns:
            list[float]: A list of volumes for the glider components.
        """
        volumetric_structure_fraction_fuselage = 0.07
        volumetric_structure_fraction_lifting_surface = 0.14
        volumetric_structure_fraction_vertical_tail = 0.27
        
        return [
            volumetric_structure_fraction_lifting_surface * self.right_wing_volume,
            volumetric_structure_fraction_lifting_surface * self.left_wing_volume,
            volumetric_structure_fraction_vertical_tail * self.vertical_tail_volume,
            volumetric_structure_fraction_lifting_surface * self.right_hor_tail_volume,
            volumetric_structure_fraction_lifting_surface * self.left_hor_tail_volume,
            volumetric_structure_fraction_fuselage * self.fuselage_volume
        ]

    @Attribute
    def list_of_x_cog(self):
        """
        Returns a list of the x-coordinate of the CG positions for various components.
        
        Returns:
            list[float]: A list of x-CG positions for each component.
        """
        return [
            self.right_wing_cog[0], self.left_wing_cog[0], self.vertical_tail_cog[0],
            self.right_hor_tail_cog[0], self.left_hor_tail_cog[0], self.fuselage_cog[0]
        ]

    @Attribute
    def list_of_masses(self):
        """
        Returns a list of the masses of each component of the glider.
        
        Returns:
            list[float]: A list of component masses (kg).
        """
        return [item * self.material_density for item in self.list_of_volumes]

    @Attribute
    def get_total_empty_mass(self):
        """
        Returns the total empty mass of the glider, including the structure and all components.
        
        Returns:
            float: The total empty mass of the glider (kg).
        """
        return sum(self.list_of_masses)

    @Attribute
    def get_moments(self):
        """
        Returns a list of the moments for each component, based on its mass and CG position.
        
        Returns:
            list[float]: A list of moments (Nm) for each component.
        """
        # Note: Datum is the nose of the fuselage (0, 0, 0)
        moment_list = []
        for i in range(len(self.list_of_masses)):
            moment_list.append(self.list_of_masses[i] * self.list_of_x_cog[i])
        return moment_list

    @Attribute
    def get_current_cog(self):
        """
        Returns the current center of gravity (CG) position of the glider, including the pilot.
        
        Returns:
            float: The current CG position relative to the MAC (Mean Aerodynamic Chord).
        """
        pilot_mass = self.current_pilot_mass
        pilot_arm = self.pilot_cog[0]
        pilot_moment = pilot_mass * pilot_arm

        total_mass = sum(self.list_of_masses) + pilot_mass
        total_moment = sum(self.get_moments) + pilot_moment
        return total_moment / total_mass

    @Attribute
    def get_fwd_limit_cog(self):
        """
        Returns the forward limit of the CG position, considering the maximum pilot mass.
        
        Returns:
            float: The forward limit of the CG position relative to the MAC.
        """
        pilot_mass_max = self.max_pilot_mass
        pilot_arm = self.pilot_cog[0]
        pilot_moment_max = pilot_mass_max * pilot_arm

        total_mass = sum(self.list_of_masses) + pilot_mass_max
        total_moment_max = sum(self.get_moments) + pilot_moment_max
        return total_moment_max / total_mass

    @Attribute
    def get_bwd_limit_cog(self):
        """
        Returns the backward limit of the CG position, considering the minimum pilot mass.
        
        Returns:
            float: The backward limit of the CG position relative to the MAC.
        """
        pilot_mass_min = self.min_pilot_mass
        pilot_arm = self.pilot_cog[0]
        pilot_moment_min = pilot_mass_min * pilot_arm

        total_mass = sum(self.list_of_masses) + pilot_mass_min
        total_moment_min = sum(self.get_moments) + pilot_moment_min
        return total_moment_min / total_mass