# Python native imports

# Python third party imports
import numpy as np
import matplotlib.pyplot as plt

# ParaPy imports
from parapy.core import Base, Input, Attribute

class ScissorPlot(Base):
    """
    A class to create a scissor plot for analyzing the stability and controllability limits
    of an aircraft based on various inputs such as stability margin, center of gravity, and 
    aerodynamic parameters.
    
    Attributes:
        SM (float): Stability margin (in percent of MAC).
        actual_x_cog (float): Actual position of the center of gravity (CG) in absolute terms.
        fwd_limit_x_cog (float): Forward limit of CG location, corresponding to maximum pilot weight.
        bwd_limit_x_cog (float): Backward limit of CG location, corresponding to minimum pilot weight.
        x_LEMAC (float): Location of the leading edge of the mean aerodynamic chord (MAC).
        x_ac (float): Location of the aerodynamic center (AC) relative to the MAC.
        s (float): Wing surface area in square meters.
        Sh_S (float): Horizontal tail surface area relative to wing surface area.
        l_h (float): Tail length (distance between AC of wing and AC of horizontal tail).
        mac (float): Mean aerodynamic chord length of the wing in meters.
        cm_ac (float): Pitching moment coefficient about the aerodynamic center.
        cl_a_min_h (float): Lift coefficient of the airplane minus the horizontal tail.
        cl_h (float): Lift coefficient of the horizontal tail.
        cl_alpha_h (float): Lift curve slope of the horizontal tail.
        cl_alpha_a_min_h (float): Lift curve slope of the airplane minus the horizontal tail.
        Vh_V (float): Ratio of velocity at the horizontal tail to velocity at the main wing.
    """
    
    # Designer input for stability margin
    SM = Input()                            # Stability margin [-]

    # Center of gravity inputs (CG)
    actual_x_cog = Input()
    fwd_limit_x_cog = Input()               # Maximum allowed pilot mass
    bwd_limit_x_cog = Input()               # Minimum allowed pilot mass

    # Leading edge of the mean aerodynamic chord location
    x_LEMAC = Input()

    x_ac = Input(0.25)                      # Location of the aerodynamic center w.r.t. MAC [-]
    s = Input()                             # Wing surface area [m^2]
    Sh_S = Input()                          # Horizontal tail surface area relative to wing surface area

    l_h = Input()                           # Tail length (distance from AC wing to AC horizontal tail) [m]
    mac = Input()                           # Mean aerodynamic chord length [m]

    # Aerodynamic parameters (from aerodynamic analysis)
    cm_ac = Input()                         # Pitching moment coefficient around AC [-]
    cl_a_min_h = Input()                    # Lift coefficient of airplane minus horizontal tail [-]
    cl_h = Input()                          # Lift coefficient of horizontal tail [-]
    cl_alpha_h = Input()                    # Lift curve slope coefficient of horizontal tail [/rad]
    cl_alpha_a_min_h = Input()              # Lift curve slope coefficient of airplane minus horizontal tail (assumed equal to cl_alpha_w) [/rad]
    Vh_V = Input(1)                         # Ratio of velocity at horizontal tail to velocity at main wing

    def convert_cog_abs_to_rel_mac(self, x_cog_abs):
        """
        Converts the absolute center of gravity (CG) position to a relative position 
        with respect to the mean aerodynamic chord (MAC).
        
        Args:
            x_cog_abs (float): The absolute center of gravity (CG) position.
        
        Returns:
            float: The relative CG position with respect to MAC.
        """
        return (x_cog_abs - self.x_LEMAC) / self.mac

    # Stability and controllability limit curves
    def x_cg_stability_limit(self, sh_s):
        """
        Calculates the stability limit for the CG position based on the horizontal tail 
        surface area ratio and other aerodynamic parameters.
        
        Args:
            sh_s (float): Horizontal tail surface area relative to wing surface area.
        
        Returns:
            float: The CG position that satisfies the stability limit.
        """
        return self.x_ac + ((self.cl_alpha_h / self.cl_alpha_a_min_h) * sh_s * (self.l_h / self.mac) * self.Vh_V ** 2) - self.SM

    def x_cg_stability_limit_no_sm(self, sh_s):
        """
        Calculates the stability limit for the CG position without considering the stability margin.
        
        Args:
            sh_s (float): Horizontal tail surface area relative to wing surface area.
        
        Returns:
            float: The CG position that satisfies the stability limit without the stability margin.
        """
        return self.x_ac + ((self.cl_alpha_h / self.cl_alpha_a_min_h) * sh_s * (self.l_h / self.mac) * self.Vh_V ** 2)

    def x_cg_controllability_limit(self, sh_s):
        """
        Calculates the controllability limit for the CG position based on the horizontal tail 
        surface area ratio and other aerodynamic parameters.
        
        Args:
            sh_s (float): Horizontal tail surface area relative to wing surface area.
        
        Returns:
            float: The CG position that satisfies the controllability limit.
        """
        return self.x_ac - (self.cm_ac / self.cl_a_min_h) + ((self.cl_h / self.cl_a_min_h) * sh_s * (self.l_h / self.mac) * self.Vh_V**2)

    # Define range for horizontal tail surface area ratio (Sh/S)
    s_h_s_range = np.linspace(0.005, 0.5, 100)

    def plot_scissor_plot(self):
        """
        Plots the scissor plot showing the stability and controllability limits, 
        along with the current CG position and limits for forward and backward CG.
        """
        plt.figure(figsize=(10, 6))
        
        # Plot stability and controllability limits
        plt.plot(self.x_cg_stability_limit(self.s_h_s_range), self.s_h_s_range, label='Stability Limit', color='blue')
        plt.plot(self.x_cg_stability_limit_no_sm(self.s_h_s_range), self.s_h_s_range, label='Stability Limit (no margin)', linestyle='--', color='blue')
        plt.plot(self.x_cg_controllability_limit(self.s_h_s_range), self.s_h_s_range, label='Controllability Limit', color='red')

        # Plot current and limit CG positions
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.actual_x_cog), color='black', linestyle='-', label=f'Current Xcg')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.fwd_limit_x_cog), color='red', linestyle='dashdot', label=f'Forward Limit Xcg (max. pilot weight)')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.bwd_limit_x_cog), color='red', linestyle='dashdot', label=f'Backward Limit Xcg (min. pilot weight)')

        # Plot current horizontal tail surface ratio
        plt.axhline(self.Sh_S, color='green', linestyle='--', label=f'Current Sh/S = {self.Sh_S:.2f}')
        
        # Set plot labels and title
        plt.xlabel('Xcg / MAC [-]')
        plt.ylabel('Sh/S [-]')
        plt.title('Scissor Plot')
        plt.grid(True)
        plt.legend()
        plt.show()