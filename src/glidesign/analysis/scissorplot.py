# Python native imports

# Python third party imports
import numpy as np
import matplotlib.pyplot as plt

# ParaPy imports
from parapy.core import Base, Input, Attribute
class ScissorPlot(Base):
    #Designer input
    SM = Input()                            #Stability margin [-]

    #Center of gravity inputs
    actual_x_cog = Input()
    fwd_limit_x_cog = Input()               #maximum allowed pilot mass
    bwd_limit_x_cog = Input()               #minimum allowed pilot mass

    # Leading edge of mean aerodynamic chord location
    x_LEMAC = Input()

    x_ac = Input(0.25)                      # Location aerodynamic centre w.r.t. MAC [-]
    s = Input()                             # Wing surface area  [m^2]
    Sh_S = Input()

    l_h = Input()                           # Tail length (distance AC wing - AC horizontal tail) [m]
    mac = Input()                           # Mean aerodynamic wing chord length [m]

    #Aerodynamic parameters (from aerodynamic analysis)
    cm_ac = Input()                         # Pitching moment coefficient around AC [-]
    cl_a_min_h = Input()                    # Lift coefficient of airplane MINUS horizontal tail [-]
    cl_h = Input()                          # Lift coefficient of horizontal tail [-]
    cl_alpha_h = Input()                    # Lift curve slope coefficient of horizontal tail [/rad]
    cl_alpha_a_min_h = Input()              # Lift curve slope coefficient of airplane MINUS horizontal tail (assumed equal to cl_alpha_w) [/rad]
    Vh_V = Input(1)                         # Ratio of velocity at the horizontal tail to velocity at the main wing

    def convert_cog_abs_to_rel_mac(self, x_cog_abs):
        #make leading edge (should technically be LE of MAC(adjust in glider.py)) the datum
        return (x_cog_abs - self.x_LEMAC)/self.mac

    #Stability and Controllability limit curves
    def x_cg_stability_limit(self, sh_s):
        return self.x_ac + ((self.cl_alpha_h / self.cl_alpha_a_min_h) * sh_s * (self.l_h / self.mac) * self.Vh_V ** 2) - self.SM

    def x_cg_stability_limit_no_sm(self, sh_s):
        return self.x_ac + ((self.cl_alpha_h / self.cl_alpha_a_min_h) * sh_s * (self.l_h / self.mac) * self.Vh_V ** 2)

    def x_cg_controllability_limit(self, sh_s):
        return self.x_ac - (self.cm_ac / self.cl_a_min_h) + ((self.cl_h / self.cl_a_min_h) * sh_s * (self.l_h /self.mac) * self.Vh_V**2)
    
    #Define Sh/S range
    s_h_s_range = np.linspace(0.005, 0.5, 100)

    #Plot the scissor plot
    def plot_scissor_plot(self):
        plt.figure(figsize=(10,6))
        plt.plot(self.x_cg_stability_limit(self.s_h_s_range), self.s_h_s_range, label='Stability Limit', color='blue')
        plt.plot(self.x_cg_stability_limit_no_sm(self.s_h_s_range), self.s_h_s_range, label='Stability Limit (no margin)', linestyle = '--', color='blue')
        plt.plot(self.x_cg_controllability_limit(self.s_h_s_range), self.s_h_s_range, label='Controllability Limit', color='red')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.actual_x_cog), color='black', linestyle='-', label=f'Current Xcg')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.fwd_limit_x_cog), color='red', linestyle='dashdot', label=f'Forward Limit Xcg (max. pilot weight)')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.bwd_limit_x_cog), color='red', linestyle='dashdot', label=f'Backward Limit Xcg (min. pilot weight)')
        plt.axhline(self.Sh_S, color='green', linestyle='--', label=f'Current Sh/S = {self.Sh_S:.2f}')
        plt.xlabel('Xcg / MAC [-]')
        plt.ylabel('Sh/S [-]')
        plt.title('Scissor Plot')
        plt.grid(True)
        plt.legend()
        plt.show()