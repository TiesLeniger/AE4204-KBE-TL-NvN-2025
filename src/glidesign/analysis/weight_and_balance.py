# Python native imports

# Python third party imports

# ParaPy imports
from parapy.core import Base, Input, Attribute
from parapy.core.validate import Range
from parapy.core.widgets import Dropdown

# Custom imports

class WeightAndBalance(Base):

    min_pilot_mass = Input()                 #Minimum allowed mass for pilot
    max_pilot_mass = Input()               #Maximum allowed pilot mass
    current_pilot_mass = Input()

    glider_structure_material = Input()

    right_wing_cog = Input()
    left_wing_cog = Input()
    vertical_tail_cog = Input()
    right_hor_tail_cog = Input()
    left_hor_tail_cog = Input()
    fuselage_cog = Input()

    right_wing_volume = Input()
    left_wing_volume = Input()
    vertical_tail_volume = Input()
    right_hor_tail_volume = Input()
    left_hor_tail_volume = Input()
    fuselage_volume = Input()

    pilot_cog = Input() #Based on center of glider CANOPY (ellipse)

    @Attribute
    def material_density(self):
        if self.glider_structure_material == "Carbon fibre":
            return 1700 #kg/m3
        elif self.glider_structure_material == "Glass fibre":
            return 2000 #kg/m3

    @Attribute
    def list_of_volumes(self):
        #Give the structure a volume that resembles the thickness of the structure
        volumetric_structure_fraction_fuselage = 0.07
        volumetric_structure_fraction_lifting_surface = 0.14
        volumetric_structure_fraction_vertical_tail = 0.27
        return [volumetric_structure_fraction_lifting_surface*self.right_wing_volume,
                volumetric_structure_fraction_lifting_surface*self.left_wing_volume,
                volumetric_structure_fraction_vertical_tail*self.vertical_tail_volume,
                volumetric_structure_fraction_lifting_surface*self.right_hor_tail_volume,
                volumetric_structure_fraction_lifting_surface*self.left_hor_tail_volume,
                volumetric_structure_fraction_fuselage*self.fuselage_volume]

    @Attribute
    def list_of_x_cog(self):
        return [self.right_wing_cog[0], self.left_wing_cog[0], self.vertical_tail_cog[0], self.right_hor_tail_cog[0], self.left_hor_tail_cog[0], self.fuselage_cog[0]]

    @Attribute
    def list_of_masses(self):
        return [item * self.material_density for item in self.list_of_volumes]

    @Attribute
    def get_total_empty_mass(self):
        return sum(self.list_of_masses)

    @Attribute
    def get_moments(self):
        #Note: datum is (0, 0, 0) (Nose of fuselage)
        moment_list = []
        for i in range(len(self.list_of_masses)):
            moment_list.append(self.list_of_masses[i] * self.list_of_x_cog[i])
        return moment_list

    @Attribute
    def get_current_cog(self):
        pilot_mass = self.current_pilot_mass
        pilot_arm = self.pilot_cog[0]
        pilot_moment = pilot_mass * pilot_arm

        total_mass = sum(self.list_of_masses) + pilot_mass
        total_moment = sum(self.get_moments) + pilot_moment
        return total_moment / total_mass

    @Attribute
    def get_fwd_limit_cog(self):
        pilot_mass_max = self.max_pilot_mass
        pilot_arm = self.pilot_cog[0]
        pilot_moment_max = pilot_mass_max * pilot_arm

        total_mass = sum(self.list_of_masses) + pilot_mass_max
        total_moment_max = sum(self.get_moments) + pilot_moment_max
        return total_moment_max / total_mass

    @Attribute
    def get_bwd_limit_cog(self):
        pilot_mass_min = self.min_pilot_mass
        pilot_arm = self.pilot_cog[0]
        pilot_moment_min = pilot_mass_min * pilot_arm

        total_mass = sum(self.list_of_masses) + pilot_mass_min
        total_moment_min = sum(self.get_moments) + pilot_moment_min
        return total_moment_min / total_mass